"""Waveforms with an explicit geocentric reference time of zero."""
import numpy as np
from pycbc.detector import Detector
from pycbc.types import FrequencySeries, TimeSeries
from pycbc.waveform import get_fd_waveform, get_td_waveform

def response(ifo, config):
    det = Detector(ifo)
    fp, fc = det.antenna_pattern(config.right_ascension, config.declination,
                                config.polarisation, config.gps)
    delay = det.time_delay_from_earth_center(config.right_ascension,
                                           config.declination, config.gps)
    return float(fp), float(fc), float(delay)

def project_td(hp, hc, ifo, config):
    fp, fc, delay = response(ifo, config)
    h = fp * hp + fc * hc
    h.start_time = float(h.start_time) + delay
    return h

def seobnrv4_template(ifo, config):
    hp, hc = get_td_waveform(approximant="SEOBNRv4_opt",
        mass1=config.detector_mass_1, mass2=config.detector_mass_2,
        spin1z=0, spin2z=0, distance=config.distance_mpc,
        inclination=config.inclination, coa_phase=0, delta_t=config.delta_t, f_lower=15)
    return project_td(hp.taper_timeseries(location="start"),
                      hc.taper_timeseries(location="start"), ifo, config)

def sine_gaussian_template(ifo, config):
    duration = 4.0
    t = (np.arange(int(duration*config.sample_rate)) / config.sample_rate) - duration/2
    # Convention: h(t)=A exp[-t^2/(2 sigma_t^2)] cos(2 pi f0 t), Q=2 pi f0 sigma_t.
    sigma_t = config.burst_q/(2*np.pi*config.burst_frequency)
    envelope = 1e-21 * np.exp(-0.5*(t/sigma_t)**2)
    hp = TimeSeries(envelope*np.cos(2*np.pi*config.burst_frequency*t),
                    delta_t=config.delta_t, epoch=-duration/2)
    hc = TimeSeries(config.burst_ellipticity*envelope*np.sin(2*np.pi*config.burst_frequency*t),
                    delta_t=config.delta_t, epoch=-duration/2)
    return project_td(hp, hc, ifo, config)

def to_reference_fd(template, data_length, delta_f):
    """FFT at physical time zero, including a TD template's relative epoch.

    Appending zeros alone discards the reference-time convention. Multiplying
    by exp(-2 pi i f t_start) restores it. The epoch is then set to zero because
    sample times are handled by the data series in matched filtering.
    """
    if isinstance(template, FrequencySeries):
        result = template.copy()
        result.resize(data_length//2 + 1)
    else:
        start = float(template.start_time)
        output = template.copy()
        if len(output) > data_length:
            raise ValueError("Template longer than analysis interval")
        output.resize(data_length)
        result = output.to_frequencyseries(delta_f=delta_f)
        result *= np.exp(-2j*np.pi*np.asarray(result.sample_frequencies)*start)
    result.start_time = 0
    return result

def imrphenomxphm_template(ifo, delta_f, length, config):
    hp, hc = get_fd_waveform(approximant="IMRPhenomXPHM",
        mass1=config.detector_mass_1, mass2=config.detector_mass_2,
        spin1x=0, spin1y=0, spin1z=0, spin2x=0, spin2y=0, spin2z=0,
        distance=config.distance_mpc, inclination=config.inclination,
        coa_phase=0, delta_f=delta_f, f_lower=15, f_ref=20,
        f_final=config.sample_rate/2)
    fp, fc, delay = response(ifo, config)
    result = fp*hp + fc*hc
    result.resize(length)
    result *= np.exp(-2j*np.pi*np.asarray(result.sample_frequencies)*delay)
    result.start_time = 0
    return result

def build_models(ifo, data_length, delta_f, config):
    return {
        "Sine Gaussian (fixed)": to_reference_fd(sine_gaussian_template(ifo,config),data_length,delta_f),
        "SEOBNRv4_opt": to_reference_fd(seobnrv4_template(ifo,config),data_length,delta_f),
        "IMRPhenomXPHM": imrphenomxphm_template(ifo,delta_f,data_length//2+1,config)}
