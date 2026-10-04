"""Official GWOSC retrieval, file integrity and data quality accounting."""
import hashlib
import json
from urllib.request import urlretrieve
import h5py
import numpy as np
from gwpy.timeseries import TimeSeries
from .config import AnalysisConfig

GWOSC_BASE = "https://gwosc.org/eventapi/json/GWTC-2.1-confident/GW190521/v4"
PREFIX = {"H1": "H-H1", "L1": "L-L1", "V1": "V-V1"}

def gwosc_url(ifo, duration=4096):
    if ifo not in PREFIX or duration != 4096:
        raise ValueError("This off-source analysis requires H1/L1/V1 4096 s files")
    return f"{GWOSC_BASE}/{PREFIX[ifo]}_GWOSC_4KHZ_R1-1242440920-4096.hdf5"

def local_path(data_dir, ifo, duration=4096):
    return data_dir / f"{ifo}_GWOSC_4KHZ_R1-1242440920-{duration}.hdf5"

def download_strain(config, duration=4096):
    config.data_dir.mkdir(parents=True, exist_ok=True)
    paths = {}
    for ifo in config.ifos:
        url = gwosc_url(ifo, duration)
        path = local_path(config.data_dir, ifo, duration)
        if not path.exists():
            partial = path.with_suffix(".part")
            print(f"Downloading {ifo} public GWOSC strain", flush=True)
            urlretrieve(url, partial)
            partial.replace(path)
        paths[ifo] = path
    return paths

def _quality(f, group, field, names_field, start, end):
    ds = f[f"quality/{group}/{field}"]
    times = ds.attrs["Xstart"] + np.arange(len(ds)) * ds.attrs["Xspacing"]
    mask = (times < end) & (times + ds.attrs["Xspacing"] > start)
    values = ds[:][mask]
    names = [x.decode() for x in f[f"quality/{group}/{names_field}"][:]]
    return {name: {"pass_seconds": int(np.count_nonzero(values & (1 << i))),
                   "total_seconds": int(len(values))}
            for i, name in enumerate(names)}

def data_manifest(config, paths):
    result = {}
    for ifo, path in paths.items():
        with path.open("rb") as stream:
            digest = hashlib.file_digest(stream, "sha256").hexdigest()
        with h5py.File(path) as f:
            intervals = {}
            for name, start, end in [
                ("psd", config.psd_start, config.psd_end),
                ("event", config.event_start, config.event_start + config.event_duration)]:
                dq = _quality(f, "simple", "DQmask", "DQShortnames", start, end)
                injections = _quality(f, "injections", "Injmask", "InjShortnames", start, end)
                intervals[name] = {"start_gps": start, "end_gps": end,
                                   "data_quality": dq, "hardware_injections": injections}
                for flag in ["DATA", "CBC_CAT1", "CBC_CAT2", "BURST_CAT1", "BURST_CAT2"]:
                    if name == "event" and flag in dq and dq[flag]["pass_seconds"] != dq[flag]["total_seconds"]:
                        raise ValueError(f"{ifo} {name} fails {flag}; inspect before analysis")
                for flag in ["NO_CBC_HW_INJ", "NO_BURST_HW_INJ"]:
                    if name == "event" and flag in injections and injections[flag]["pass_seconds"] != injections[flag]["total_seconds"]:
                        raise ValueError(f"{ifo} {name} contains a hardware injection")
            result[ifo] = {"url": gwosc_url(ifo), "bytes": path.stat().st_size,
                           "sha256": digest, "intervals": intervals}
    return result

def load_segments(config, paths):
    long_data, event_data = {}, {}
    for ifo, path in paths.items():
        # Resample an expanded interval, then discard the filter boundaries.
        series = TimeSeries.read(path, format="hdf5.gwosc",
            start=config.psd_start - 8,
            end=config.event_start + config.event_duration + 8).resample(config.sample_rate)
        if not np.isfinite(series.value).all():
            raise ValueError(f"Non-finite strain for {ifo}")
        long_data[ifo] = series.crop(config.psd_start, config.psd_end)
        event_data[ifo] = series.crop(config.event_start, config.event_start + config.event_duration)
    return long_data, event_data


def valid_psd_windows(config, paths):
    """Keep 8 s windows passing CBC/burst CAT2 flags with 1 s veto padding."""
    output = {}
    stride = config.psd_fftlength - config.psd_overlap
    for ifo,path in paths.items():
        accepted = []
        with h5py.File(path) as f:
            for start in range(config.psd_start, config.psd_end-config.psd_fftlength+1, stride):
                end = start + config.psd_fftlength
                dq = _quality(f,"simple","DQmask","DQShortnames",start-1,end+1)
                inj = _quality(f,"injections","Injmask","InjShortnames",start-1,end+1)
                required = [dq[n] for n in ["DATA","CBC_CAT1","CBC_CAT2","BURST_CAT1","BURST_CAT2"] if n in dq]
                required += [inj[n] for n in ["NO_CBC_HW_INJ","NO_BURST_HW_INJ"] if n in inj]
                if all(v["pass_seconds"]==v["total_seconds"] and v["total_seconds"]>0 for v in required):
                    accepted.append(start)
        if len(accepted)<16:
            raise ValueError(f"Insufficient good PSD windows for {ifo}")
        output[ifo]=accepted
    return output
