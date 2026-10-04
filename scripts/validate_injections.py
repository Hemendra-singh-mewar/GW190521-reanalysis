#!/usr/bin/env python3
"""Conditional amplitude recovery in 128 synthetic Gaussian-noise realisations.

This validates the known-shape, known-time likelihood normalisation. It does
not assess astrophysical parameter bias or posterior coverage.
"""
from dataclasses import replace
import json
import numpy as np
from pycbc.filter import matched_filter, sigma
from pycbc.psd import interpolate
from pycbc.types import FrequencySeries
from gw190521.config import AnalysisConfig
from gw190521.data import download_strain,load_segments,valid_psd_windows
from gw190521.conditioning import estimate_psds,analysis_strain
from gw190521.waveforms import build_models,sine_gaussian_template,to_reference_fd
from gw190521.filtering import filter_series,network_summary

c=AnalysisConfig()
paths=download_strain(c)
long,event=load_segments(c,paths)
psds=estimate_psds(long,c,valid_psd_windows(c,paths))
n=c.event_duration*c.sample_rate
nf=n//2+1;df=1/c.event_duration
f=np.arange(nf)*df
injection_offset=16.375
injection_time=c.event_start+injection_offset
trial_config=replace(c,gps=injection_time)
models={ifo:build_models(ifo,n,df,c) for ifo in c.ifos}
rng=np.random.default_rng(20261002)
output={"seed":20261002,"realisations_per_model":128,"network_optimal_snr":20,
        "injection_time_gps":injection_time,"scope":"Known template shape, geometry and time; real common amplitude only; ideal Gaussian noise from measured PSDs.","models":{}}
for model in ['Sine Gaussian (fixed)','SEOBNRv4_opt','IMRPhenomXPHM']:
    templates={};norms={};spectra={}
    for ifo in c.ifos:
        p=interpolate(psds[ifo],df)
        h=models[ifo][model]
        norms[ifo]=float(sigma(h,psd=p,low_frequency_cutoff=c.f_low,high_frequency_cutoff=c.f_high))
        spectra[ifo]=p
    scale=20/np.linalg.norm(list(norms.values()))
    for ifo in c.ifos:
        templates[ifo]=models[ifo][model]*scale
        norms[ifo]*=scale
    residuals=[]
    for trial in range(128):
        numerator=0
        for ifo in c.ifos:
            p=spectra[ifo]
            noise=(rng.normal(size=nf)+1j*rng.normal(size=nf))*np.sqrt(np.asarray(p)/(4*df))
            noise[[0,-1]]=0
            signal=np.asarray(templates[ifo])*np.exp(-2j*np.pi*f*injection_offset)
            d=FrequencySeries(noise+signal,delta_f=df).to_timeseries(delta_t=c.delta_t)
            d.start_time=c.event_start
            z=matched_filter(templates[ifo],d,psd=p,low_frequency_cutoff=c.f_low,
                             high_frequency_cutoff=c.f_high)
            numerator += norms[ifo]*float(np.real(z[int(round(injection_offset*c.sample_rate))]))
        amplitude=numerator/(20**2)
        residuals.append((amplitude-1)*20)
    res=np.asarray(residuals)
    output['models'][model]={
        'mean_standardised_amplitude_error':float(res.mean()),
        'sd_standardised_amplitude_error':float(res.std(ddof=1)),
        'fraction_within_1_sigma':float(np.mean(np.abs(res)<=1)),
        'fraction_within_2_sigma':float(np.mean(np.abs(res)<=2))}
# Isolate the epoch bug using exactly the same new sine Gaussian shape and PSD.
correct_z={};legacy_z={};norms={}
for ifo in c.ifos:
    d=analysis_strain(event[ifo]);h=sine_gaussian_template(ifo,c)
    correct_z[ifo],norms[ifo]=filter_series(h,d,psds[ifo],c)
    wrong=h.copy();wrong.resize(len(d));wrong=wrong.to_frequencyseries();wrong.start_time=0
    legacy_z[ifo],_=filter_series(wrong,d,psds[ifo],c)
output['epoch_ablation']={
    'correct_reference':network_summary(correct_z,norms,c)[0],
    'epoch_discarded':network_summary(legacy_z,norms,c)[0]}
c.results_dir.mkdir(exist_ok=True)
(c.results_dir/'injection_validation.json').write_text(json.dumps(output,indent=2)+'\n')
print(json.dumps(output,indent=2))
