#!/usr/bin/env python3
"""Run the public GW190521 model-comparison analysis."""

from gw190521.pipeline import run_analysis


if __name__ == "__main__":
    snr, matches = run_analysis()
    print("\nMatched-filter SNR summary\n")
    print(snr.round(3).to_string())
    print("\nPSD weighted model matches (H1 weighting)\n")
    print(matches.round(3).to_string())
