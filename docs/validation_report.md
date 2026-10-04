# Validation report

Date: 2 October 2026. Software: version 0.2.0.

## Completed checks

- End-to-end analysis on the official three 4096 s GWOSC records.
- Full-event DATA, CBC CAT1/CAT2, BURST CAT1/CAT2 and relevant hardware-injection flags.
- Off-source PSD screening: 124 H1, 108 L1 and 127 V1 accepted windows of 127 candidates.
- Fresh-kernel execution of all seven notebook code cells, with zero errors.
- Eleven passing regression cases, including analytic burst injections, independent
  time placement of CBC injections, reference-time invariance under padding,
  frequency-domain detector delay, network-statistic ordering and the Gaussian likelihood.
- 128 conditional amplitude trials for each of three fixed templates, 384 in total.
- Frequency-band and PSD-duration checks with template parameters held fixed.

## Results

| Template | Independent peak bound | Common time | Fixed response coherent |
|---|---:|---:|---:|
| Sine Gaussian (fixed) | 14.314 | 14.175 | 6.862 |
| SEOBNRv4_opt | 14.876 | 14.481 | 8.487 |
| IMRPhenomXPHM | 14.170 | 13.949 | 8.301 |
| Sine Gaussian (grid) | 14.640 | 14.466 | 7.743 |

`results/injection_validation.json` contains the seed, tested assumptions and
all simulation summary values. `results/analysis_metadata.json` contains the
input URLs, hashes, quality counts, accepted window start times, detector
responses, configuration and package versions. The test suite requires no raw
strain or network connection.

## Corrected development errors

The previous reconstruction discarded the relative epoch of padded time-domain
templates. The burst was consequently evaluated at the wrong lag. It also called
the quadrature sum of independently selected detector peaks a coincident network
statistic without enforcing common arrival-time consistency. Both errors are
corrected and directly tested. The burst Q convention is now explicit, and the
PSD used for plots is the same off-source PSD used for matched filtering.

## What these checks do not establish

No mass, spin, distance or sky posterior has been inferred. The conditional
amplitude simulation assumes known template shape, phase, geometry and time in
ideal Gaussian noise. It is not an astrophysical injection campaign or a
convergence study. Full calibration/PSD uncertainty, coherent glitch treatment,
noise background estimation and Bayesian model selection are absent. The
fixed-response statistic uses an illustrative, unfitted geometry and should
not be used to reject the coherence of the real event.

The public note and README state these limitations explicitly.
