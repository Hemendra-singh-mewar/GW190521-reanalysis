"""Matched filters and explicitly defined network summary statistics."""
import numpy as np
from pycbc.filter import match, matched_filter, sigma
from pycbc.psd import interpolate
from .waveforms import to_reference_fd

def filter_series(template, data, psd, config):
    h = to_reference_fd(template, len(data), data.delta_f)
    p = interpolate(psd, data.delta_f)
    args = dict(psd=p, low_frequency_cutoff=config.f_low,
                high_frequency_cutoff=config.f_high)
    z = matched_filter(h,data,**args).crop(4,4)
    return z, float(sigma(h,**args))

def network_summary(series, norms, config):
    reference = next(iter(series.values()))
    time = np.asarray(reference.sample_times)
    selection = np.abs(time-config.gps) <= config.on_source_half_width
    if not selection.any():
        raise ValueError("Event window outside analysis interval")
    time = time[selection]
    z = np.asarray([np.asarray(series[ifo])[selection] for ifo in config.ifos])
    weights = np.asarray([norms[ifo] for ifo in config.ifos])
    peaks = np.max(np.abs(z),axis=1)
    peak_idx = np.argmax(np.abs(z),axis=1)
    q = np.sqrt(np.sum(np.abs(z)**2,axis=0))
    # Common complex coefficient, fixed detector amplitudes, phases and delays.
    coherent = np.abs(np.sum(weights[:,None]*z,axis=0)) / np.linalg.norm(weights)
    iq, ic = int(np.argmax(q)), int(np.argmax(coherent))
    row = dict(zip(config.ifos,peaks.astype(float)))
    row.update(Independent_peak_bound=float(np.linalg.norm(peaks)),
               Time_aligned=float(q[iq]), Fixed_response_coherent=float(coherent[ic]))
    timing = {f"{ifo}_peak_geocentre_offset_s":float(time[i]-config.gps)
              for ifo,i in zip(config.ifos,peak_idx)}
    timing.update(time_aligned_offset_s=float(time[iq]-config.gps),
                  coherent_offset_s=float(time[ic]-config.gps))
    curve = {"offset_s":time-config.gps,"time_aligned":q,"coherent":coherent}
    return row, timing, curve

def model_match(first, second, data_length, delta_f, psd, config):
    a = to_reference_fd(first,data_length,delta_f)
    b = to_reference_fd(second,data_length,delta_f)
    value,_ = match(a,b,psd=interpolate(psd,delta_f),
        low_frequency_cutoff=config.f_low,high_frequency_cutoff=config.f_high,
        subsample_interpolation=True)
    return float(value)

def log_likelihood_ratio(data, template, psd, config):
    """Gaussian log likelihood ratio for already aligned real strains.

    ln L(h)-ln L(0)=(d|h)-0.5(h|h). No amplitude/time maximisation here.
    """
    d = data.to_frequencyseries() if hasattr(data,"to_frequencyseries") else data
    h = template.to_frequencyseries() if hasattr(template,"to_frequencyseries") else template
    p = interpolate(psd,d.delta_f)
    f = np.asarray(d.sample_frequencies)
    mask=(f>=config.f_low)&(f<config.f_high)
    dh = 4*d.delta_f*np.real(np.sum(np.asarray(d)[mask].conj()*np.asarray(h)[mask]/np.asarray(p)[mask]))
    hh = 4*d.delta_f*np.sum(np.abs(np.asarray(h)[mask])**2/np.asarray(p)[mask])
    return float(dh-0.5*hh)
