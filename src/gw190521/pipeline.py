"""Reproduce the complete fixed morphology study and its audit trail."""
from dataclasses import asdict, replace
from importlib.metadata import version
import json
import platform
import numpy as np
import pandas as pd
from .config import AnalysisConfig
from .data import download_strain, load_segments, data_manifest, valid_psd_windows
from .conditioning import estimate_psds, conditioned_event, analysis_strain
from .waveforms import build_models, sine_gaussian_template, to_reference_fd, response
from .filtering import filter_series, network_summary, model_match
from .plotting import save_all


def evaluate(models, strain, psds, config):
    rows, times, curves = {}, {}, {}
    for model in models[config.ifos[0]]:
        zs, norms = {}, {}
        for ifo in config.ifos:
            zs[ifo], norms[ifo] = filter_series(models[ifo][model],strain[ifo],psds[ifo],config)
        rows[model], times[model], curves[model] = network_summary(zs,norms,config)
    return pd.DataFrame.from_dict(rows,orient="index"),times,curves


def scan_sine_gaussian(strain, psds, config):
    records=[]
    best=None
    for frequency in np.arange(40.,90.1,2.):
        for quality in np.arange(2.,20.1,1.):
            c=replace(config,burst_frequency=float(frequency),burst_q=float(quality))
            zs,norms={},{}
            for ifo in c.ifos:
                h=sine_gaussian_template(ifo,c)
                zs[ifo],norms[ifo]=filter_series(h,strain[ifo],psds[ifo],c)
            row,timing,_=network_summary(zs,norms,c)
            records.append(dict(frequency_hz=frequency,q=quality,**row,**timing))
            if best is None or row["Time_aligned"]>best[0]:
                best=(row["Time_aligned"],c)
    return pd.DataFrame(records),best[1]


def run_analysis(config=None, *, make_qtransforms=True, scan=True):
    config=config or AnalysisConfig()
    config.results_dir.mkdir(parents=True,exist_ok=True)
    config.figures_dir.mkdir(parents=True,exist_ok=True)
    paths=download_strain(config)
    manifest=data_manifest(config,paths)
    windows=valid_psd_windows(config,paths)
    long_data,event_data=load_segments(config,paths)
    psds=estimate_psds(long_data,config,windows)
    strain={ifo:analysis_strain(d) for ifo,d in event_data.items()}
    models={ifo:build_models(ifo,len(d),d.delta_f,config) for ifo,d in strain.items()}
    scan_table=None
    if scan:
        print("Scanning 494 declared sine Gaussian shapes",flush=True)
        scan_table,best=scan_sine_gaussian(strain,psds,config)
        scan_table.to_csv(config.results_dir/'sine_gaussian_scan.csv',index=False)
        for ifo,d in strain.items():
            models[ifo]["Sine Gaussian (grid)"]=to_reference_fd(sine_gaussian_template(ifo,best),len(d),d.delta_f)
    table,timing,curves=evaluate(models,strain,psds,config)
    table.index.name="Model"
    table.to_csv(config.results_dir/'snr_summary.csv',float_format='%.9f')
    pd.DataFrame(timing).T.to_csv(config.results_dir/'timing_summary.csv',float_format='%.9f')
    names=list(models['H1'])
    matches={}
    match_rows=[]
    for ifo,d in strain.items():
        mt=pd.DataFrame(index=names,columns=names,dtype=float)
        for first in names:
            for second in names:
                mt.loc[first,second]=model_match(models[ifo][first],models[ifo][second],len(d),d.delta_f,psds[ifo],config)
                match_rows.append(dict(detector=ifo,model_a=first,model_b=second,match=mt.loc[first,second]))
        matches[ifo]=mt
    pd.DataFrame(match_rows).to_csv(config.results_dir/'model_matches_all_detectors.csv',index=False)
    matches['H1'].to_csv(config.results_dir/'model_matches.csv',float_format='%.9f')
    sensitivity=[]
    for low,high in [(20.,128.),(25.,256.),(30.,256.)]:
        c=replace(config,f_low=low,f_high=high)
        tab,_,_=evaluate(models,strain,psds,c)
        for name,row in tab.iterrows():
            sensitivity.append(dict(setting=f'band {low:g} to {high:g} Hz',model=name,**row.to_dict()))
    # Keep the chosen templates fixed when changing the PSD estimate.
    c=replace(config,psd_duration=256)
    w=valid_psd_windows(c,paths)
    psd_short=estimate_psds(long_data,c,w)
    tab,_,_=evaluate(models,strain,psd_short,c)
    for name,row in tab.iterrows():
        sensitivity.append(dict(setting='nearest 256 s PSD',model=name,**row.to_dict()))
    pd.DataFrame(sensitivity).to_csv(config.results_dir/'sensitivity.csv',index=False)
    np.savez_compressed(config.results_dir/'filter_curves.npz',
        **{name+'__'+key:value for name,curve in curves.items() for key,value in curve.items()})
    metadata=dict(config={k:str(v) if hasattr(v,'as_posix') else v for k,v in asdict(config).items()},
        environment={p:version(p) for p in ['numpy','scipy','pycbc','lalsuite','gwpy','gwosc','astropy','h5py','matplotlib','pandas']},
        python=platform.python_version(),platform=platform.platform(),
        mass_1_detector_msun=config.detector_mass_1,mass_2_detector_msun=config.detector_mass_2,
        psd_intervals={ifo:starts for ifo,starts in windows.items()},
        psd_estimator='median of Hann periodograms, constant detrending, PyCBC median bias correction',
        data=manifest,
        responses={ifo:dict(zip(['F_plus','F_cross','delay_from_geocentre_s'],response(ifo,config))) for ifo in config.ifos},
        grid_best=dict(frequency_hz=best.burst_frequency,q=best.burst_q) if scan else None,
        scope='Fixed waveform shapes and illustrative fixed geometry; no posterior, Bayes factor, detection significance or source-parameter recovery.')
    (config.results_dir/'analysis_metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    conditioned=conditioned_event(event_data,config,psds)
    save_all(table,matches['H1'],psds,event_data,conditioned,curves,scan_table,config,make_qtransforms)
    return table,matches['H1']
