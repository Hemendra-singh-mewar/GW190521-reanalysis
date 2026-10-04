# Audit of `HS004.ipynb`

Historical source audit from 22 September 2026, retained for provenance.
The corrections made on 2 October 2026 are recorded in `validation_report.md`.

This review records how the 210-cell coursework notebook was translated into
the public project. It protects the useful scientific work while separating it
from results that cannot be reproduced from the displayed execution state.

| Original material | Decision | Public version |
|---|---|---|
| GW190521 motivation and intermediate-mass black-hole context | Rewrite | Concise introduction with primary references |
| Statement that GWs arise from a changing dipole | Remove | Leading radiation described through the mass quadrupole |
| H1, L1 and V1 GWOSC data retrieval | Keep | Official GWTC 2.1 v4 HDF5 URLs, downloaded on demand |
| 2048 Hz resampling | Keep | Correctly justified by the 20 to 256 Hz analysis band |
| Statement that the event has no frequencies above 1 Hz | Remove | Event described as a short signal at tens of hertz |
| Welch PSD | Keep and clarify | 512 s median Welch estimate; no claim that Welch removes nonstationarity |
| Q transforms | Keep | Common 20 to 256 Hz display for all detectors |
| Detector-specific bands chosen by eye | Remove | One common, explicit band |
| `d = 5300 # distance in parsecs` | Correct | Distance always named `distance_mpc` |
| `z=1.82` terminology | Correct | `z` and the factor `1+z` are kept distinct |
| Direct comparison of one polarisation with detector strain | Replace | Plus and cross projected through detector responses |
| `scipy.signal.correlate` as model overlap | Remove | Normalised PSD weighted PyCBC `match` |
| CBC matched filtering | Keep and rebuild | Standard PyCBC matched filter with an interpolated PSD |
| Per-detector `p=0` claims | Remove | No zero-probability claim from finite backgrounds |
| `GW150914` text in GW190521 output | Remove | Event identity defined once in configuration |
| Network likelihood idea | Retain for future extension | Version 0.2 adds a validated Gaussian likelihood; source-parameter inference remains future work |
| Unbounded positive-mass and distance priors | Remove | Future inference must use declared finite ranges |
| Burst cross polarisation multiplied by zero | Correct | Nonzero quadrature cross polarisation |
| Optimiser result later replaced by a manual array | Remove | No hidden or manually inserted fit vector |
| MCMC that crashes, followed by reported mass estimates | Remove | No posterior claims from failed execution |
| 100 walkers for 50 steps | Remove | Converged inference listed as future work |
| Claims that lensing/eccentricity were dismissed | Remove | Alternative interpretations attributed to the literature |

## Reproducibility finding

The original notebook has 95 code cells with strongly non-monotonic execution
counts, including cells numbered 85, 4, 107, 26 and 153 in document order.
Several cells have no execution count while retaining outputs. The public
version therefore rebuilds the analysis from source modules rather than trying
to preserve hidden kernel state.

## Material not redistributed

The coursework brief and marking instructions are not reproduced in the public
notebook. Any code derived substantially from private teaching material should
be checked and attributed before release. The reconstructed implementation here
uses public GWOSC, GWPy and PyCBC interfaces and newly written analysis code.
