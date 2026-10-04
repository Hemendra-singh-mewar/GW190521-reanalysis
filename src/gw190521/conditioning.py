"""Use the same off-source PSD for filtering and visualisation."""
import numpy as np
from scipy.signal.windows import tukey
from gwpy.frequencyseries import FrequencySeries as GwpyFrequencySeries
from pycbc.psd import welch, interpolate
from .config import AnalysisConfig

def estimate_psds(long_data, config, window_starts):
    from scipy.signal import periodogram
    from pycbc.psd import median_bias
    from pycbc.types import FrequencySeries
    result = {}
    for ifo,series in long_data.items():
        spectra=[]
        for start in window_starts[ifo]:
            segment = series.crop(start,start+config.psd_fftlength)
            f,p = periodogram(segment.value, fs=config.sample_rate, window="hann",
                              detrend="constant", scaling="density")
            spectra.append(p)
        estimate=np.median(spectra,axis=0)/median_bias(len(spectra))
        result[ifo]=FrequencySeries(estimate,delta_f=f[1]-f[0])
    return result

def gwpy_asd(psd):
    return GwpyFrequencySeries(np.sqrt(np.asarray(psd)), f0=0, df=psd.delta_f)

def conditioned_event(event_data, config, psds):
    return {ifo: series.whiten(asd=gwpy_asd(psds[ifo]), fduration=4,
                               highpass=config.f_low).bandpass(config.f_low, config.f_high)
                .crop(config.gps-0.5, config.gps+0.5)
            for ifo, series in event_data.items()}

def analysis_strain(series):
    output = series.to_pycbc()
    output -= np.mean(output.numpy())
    # One-second cosine taper at each end of the 32 s interval.
    output *= tukey(len(output), alpha=2.0 / output.duration)
    return output
