"""Regression tests targeting timing, projection and statistical normalisation."""
from dataclasses import replace
import numpy as np
import pytest
from pycbc.filter import sigma
from pycbc.types import FrequencySeries, TimeSeries
from gw190521.config import AnalysisConfig
from gw190521.waveforms import (sine_gaussian_template, seobnrv4_template,
    to_reference_fd, response, imrphenomxphm_template)
from gw190521.filtering import filter_series, network_summary, model_match, log_likelihood_ratio

@pytest.fixture
def setup():
    c=AnalysisConfig()
    n=c.event_duration*c.sample_rate
    psd=FrequencySeries(np.full(n//2+1,1e-44),delta_f=1/c.event_duration)
    return c,n,psd

def test_detector_frame_mass_conversion(setup):
    c,_,_=setup
    assert c.detector_mass_1 == pytest.approx(153.504)
    assert c.detector_mass_2 == pytest.approx(89.232)

@pytest.mark.parametrize('ifo',['H1','L1','V1'])
def test_analytic_burst_arrival_and_snr(setup,ifo):
    c,n,p=setup
    # Directly evaluate the injected detector signal on absolute physical times.
    # This does not use the template FFT or epoch conversion to place the signal.
    fp,fc,delay=response(ifo,c)
    physical=np.arange(n)*c.delta_t + (c.event_start-c.gps)-delay
    width=c.burst_q/(2*np.pi*c.burst_frequency)
    h=1e-21*np.exp(-.5*(physical/width)**2)*(fp*np.cos(2*np.pi*c.burst_frequency*physical)
                                            +fc*np.sin(2*np.pi*c.burst_frequency*physical))
    d=TimeSeries(h,delta_t=c.delta_t,epoch=c.event_start)
    z,norm=filter_series(sine_gaussian_template(ifo,c),d,p,c)
    idx=int(np.argmax(np.abs(z.numpy())))
    assert abs(float(z.sample_times[idx])-c.gps)<=c.delta_t
    assert abs(z[idx])/norm == pytest.approx(1,rel=3e-4)

@pytest.mark.parametrize('ifo',['H1','L1','V1'])
def test_cbc_arrival_from_independent_time_placement(setup,ifo):
    c,n,p=setup
    h=seobnrv4_template(ifo,c)
    physical=np.arange(n)*c.delta_t+(c.event_start-c.gps)
    y=np.interp(physical,np.asarray(h.sample_times),np.asarray(h),left=0,right=0)
    d=TimeSeries(y,delta_t=c.delta_t,epoch=c.event_start)
    z,norm=filter_series(h,d,p,c)
    idx=int(np.argmax(np.abs(z.numpy())))
    assert abs(float(z.sample_times[idx])-c.gps)<=c.delta_t
    assert abs(z[idx])/norm == pytest.approx(1,rel=.01)

def test_frequency_delay_is_physical(setup):
    c,n,p=setup
    from pycbc.waveform import get_fd_waveform
    hp,hc=get_fd_waveform(approximant='IMRPhenomXPHM',mass1=c.detector_mass_1,
        mass2=c.detector_mass_2,spin1x=0,spin1y=0,spin1z=0,spin2x=0,spin2y=0,spin2z=0,
        distance=c.distance_mpc,inclination=c.inclination,coa_phase=0,
        delta_f=p.delta_f,f_lower=15,f_ref=20,f_final=c.sample_rate/2)
    for ifo in c.ifos:
        fp,fc,delay=response(ifo,c)
        h=imrphenomxphm_template(ifo,p.delta_f,n//2+1,c)
        unshifted=fp*hp+fc*hc
        f=np.asarray(h.sample_frequencies)
        use=(f>25)&(f<200)&(np.abs(np.asarray(unshifted))>0)
        ratio=np.asarray(h)[use]/np.asarray(unshifted)[use]
        slope=np.polyfit(f[use],np.unwrap(np.angle(ratio)),1)[0]
        assert -slope/(2*np.pi)==pytest.approx(delay,abs=1e-10)

def test_reference_conversion_is_invariant_to_zero_padding(setup):
    c,n,p=setup
    h=sine_gaussian_template('H1',c)
    padded=TimeSeries(np.r_[np.zeros(3072),h.numpy(),np.zeros(1500)],
        delta_t=c.delta_t,epoch=float(h.start_time)-3072*c.delta_t)
    a=to_reference_fd(h,n,p.delta_f);b=to_reference_fd(padded,n,p.delta_f)
    assert np.max(np.abs(a.numpy()-b.numpy()))/np.max(np.abs(a.numpy()))<1e-6
    assert model_match(a,b,n,p.delta_f,p,c)==pytest.approx(1,abs=1e-8)

def test_joint_statistics_reject_separate_peak_sum(setup):
    c,n,p=setup
    x=np.arange(n)*c.delta_t+(c.event_start-c.gps)
    zs={ifo:TimeSeries(np.exp(-((x-shift)/.005)**2).astype(complex),
        delta_t=c.delta_t,epoch=c.event_start) for ifo,shift in zip(c.ifos,[-.1,0,.1])}
    row,_,_=network_summary(zs,dict.fromkeys(c.ifos,1.),c)
    assert row['Independent_peak_bound']>1.7
    assert row['Time_aligned']<1.01
    assert row['Fixed_response_coherent']<row['Time_aligned']

def test_likelihood_normalisation_and_maximum(setup):
    c,n,p=setup
    h=to_reference_fd(sine_gaussian_template('L1',c),n,p.delta_f)
    norm=sigma(h,psd=p,low_frequency_cutoff=c.f_low,high_frequency_cutoff=c.f_high)
    ll=log_likelihood_ratio(h,h,p,c)
    assert ll==pytest.approx(.5*norm**2,rel=1e-10)
    for amplitude in [.5,1.5]:
        assert log_likelihood_ratio(h,amplitude*h,p,c)<ll
