#!/usr/bin/env python3
"""Build a clean notebook whose results are generated rather than typed in."""
from pathlib import Path
import nbformat as nbf
r=Path(__file__).resolve().parents[1]
m=lambda t:nbf.v4.new_markdown_cell(t.strip())
c=lambda t:nbf.v4.new_code_cell(t.strip())
cells=[m(r'''
# Revisiting GW190521 with open data
## Waveform timing, noise weighting and network consistency

**Hemendra Singh · version 0.2.2 · revised 3 October 2026**

The work was carried out in 2024 at the Gravity Exploration Institute, School of Physics and Astronomy, Cardiff University.

This notebook compares fixed compact binary templates and sine Gaussian bursts
using public H1, L1 and V1 data. It accompanies the complete methods note in
`paper/`. Its scope is signal morphology and implementation validation.

No source-parameter posterior, Bayes factor or detection significance is claimed.
The CBC spins are zero and the sky/orientation is illustrative and unfitted.
A template family's ability to model precession is not evidence for precession
in this configuration.
'''),c('''
from pathlib import Path
import os, json
if Path.cwd().name == 'notebooks':
    os.chdir(Path.cwd().parent)
import pandas as pd
from IPython.display import display, Image
from gw190521.config import AnalysisConfig
from gw190521.pipeline import run_analysis
config = AnalysisConfig()
print('Event GPS:', config.gps)
print('Detector-frame masses:', config.detector_mass_1, config.detector_mass_2)
print('Luminosity distance [Mpc]:', config.distance_mpc)
'''),m(r'''
## Data and noise model

The first run downloads approximately 389 MB of official GWOSC 4096 s strain.
The analysis uses 32 s of event data resampled to 2048 Hz. The noise estimate uses
a separate 512 s interval ending 32.4 s before the catalogue reference time.
Eight-second Hann periodograms are retained only when their quality flags pass,
with a one-second guard around each window. The median PSD has a finite-sample
bias correction. It does not make detector noise stationary or Gaussian.

A common 20 to 256 Hz band is used. Filtering operates on tapered strain with
PSD weighting. Plots and Q transforms are whitened separately with the same PSD.
'''),c('''
snr_table, match_table = run_analysis(config)
display(snr_table.round(3))
'''),c('''
metadata = json.loads(Path('results/analysis_metadata.json').read_text())
print('Accepted PSD windows:', {k: len(v) for k, v in metadata['psd_intervals'].items()})
display(pd.DataFrame(metadata['responses']).T)
display(Image(filename='figures/detector_psd.png'))
display(Image(filename='figures/qtransform.png'))
'''),m(r'''
## Reference times and burst convention

The burst uses $h_+=A e^{-t^2/(2\sigma_t^2)}\cos(2\pi f_0t)$ and
$h_\times=\epsilon A e^{-t^2/(2\sigma_t^2)}\sin(2\pi f_0t)$, with
$Q=2\pi f_0\sigma_t$. Here $\epsilon=1$ gives circular polarisation.
The fixed shape has $f_0=60$ Hz and $Q=8$.

A time-domain template starting at relative time $t_s$ is converted as
$\tilde h_{\rm ref}(f)=\tilde h_{\rm array}(f)e^{-2\pi i f t_s}$.
Dropping this phase after padding shifts the burst's reference by about two
seconds. Both time and frequency-domain templates include the detector delay.
All filter outputs can then be compared at a common geocentric reference time.
'''),m(r'''
## Three different network quantities

For normalised complex filter output $z_I$ and template norm $\sigma_I$:

$$\rho_{\rm bound}=\sqrt{\sum_I\max_t|z_I(t)|^2},\qquad
\rho_{\rm time}=\max_t\sqrt{\sum_I|z_I(t)|^2},$$
$$\rho_{\rm fixed}=\max_t\frac{|\sum_I\sigma_I z_I(t)|}{\sqrt{\sum_I\sigma_I^2}}.$$

The first is an independent-peak upper bound, not a coincidence test. The
second fixes the common reference time but permits independent detector phases
and amplitudes. The third fits one common complex coefficient to fixed relative
responses. The latter is a restricted coherent projection, not a physical
source-parameter inference. An overall quadrature phase need not correspond to
one orbital phase for a multimode waveform.
'''),c('''
display(Image(filename='figures/snr_comparison.png'))
display(Image(filename='figures/network_timeseries.png'))
display(pd.read_csv('results/timing_summary.csv', index_col=0))
'''),m(r'''
## Waveform matches and morphology grid

The inner product is $(a|b)=4\operatorname{Re}\int\tilde a^*\tilde b/S_n\,df$.
The normalised match is maximised over time and filter quadrature phase. The
matrix below uses H1 noise weighting; all detector weightings are also saved.

The exploratory burst grid has 494 shapes: 40 to 90 Hz in 2 Hz steps and integer
Q from 2 to 20. It maximises the common-time statistic on the observed data, so
it has additional freedom and an uncalibrated trials effect. It is not a
posterior or an equal-complexity comparison with the fixed CBC points.
'''),c('''
display(match_table.round(3))
display(Image(filename='figures/match_matrix.png'))
print('Best sampled burst shape:', metadata['grid_best'])
display(Image(filename='figures/sine_gaussian_scan.png'))
'''),m('''
## Validation with known injections

The regression suite tests independent analytic burst timing, independently
placed time-domain CBC injections, frequency-domain delays, padding invariance,
network ordering and likelihood normalisation. The following script also runs
128 ideal Gaussian-noise realisations for each of three shapes. Only a real
common amplitude is estimated, at known shape, phase, time and geometry.
These checks do not assess astrophysical posterior coverage or parameter bias.
'''),c('''
import runpy, contextlib, io
with contextlib.redirect_stdout(io.StringIO()):
    runpy.run_path('scripts/validate_injections.py', run_name='__main__')
validation = json.loads(Path('results/injection_validation.json').read_text())
display(pd.DataFrame(validation['models']).T)
display(pd.DataFrame(validation['epoch_ablation']).T)
'''),m('''
## Sensitivity and interpretation

The next table changes the analysis band and the length of the PSD interval.
Template parameters, including the selected burst shape, remain fixed.
A correctly timed burst reaches a statistic similar in scale to the fixed CBC
examples. This establishes a useful morphology comparison; it does not show
that a burst description and a physical binary model have equal evidence.

The low fixed-response coherent values reflect the imposed, unfitted geometry.
They do not demonstrate incoherence of GW190521. The catalogue SNR is a scale
reference, not a pass/fail target for this different analysis.
'''),c('''
sensitivity = pd.read_csv('results/sensitivity.csv')
display(sensitivity[['setting', 'model', 'Time_aligned', 'Fixed_response_coherent']].round(3))
'''),m('''
## Reproducibility and limitations

Every number above is generated by code. Exact input addresses, hashes, quality
counts, accepted PSD windows and package versions are recorded in
`results/analysis_metadata.json`. No raw detector strain is committed.

A future astrophysical inference study still needs a physical joint likelihood
explored over source parameters, finite priors, converged sampling, injections
spanning its parameter region, uncertainty treatment and an official joint
posterior comparison. None of those tasks is implied complete here.

References and detailed derivations appear in the accompanying methods note.
Data release: https://gwosc.org/eventapi/html/GWTC-2.1-confident/GW190521/v4/

Acknowledgements: https://gwosc.org/acknowledgement/
''')]
nb=nbf.v4.new_notebook(cells=cells,metadata={'kernelspec':{'display_name':'Python 3','language':'python','name':'python3'},'language_info':{'name':'python','version':'3.12'}})
nbf.write(nb,r/'notebooks/GW190521_reanalysis.ipynb')
print('Notebook generated')
