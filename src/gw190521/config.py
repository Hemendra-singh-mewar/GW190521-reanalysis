"""Explicit settings for the public fixed morphology study."""
from dataclasses import dataclass, field
from pathlib import Path

@dataclass(frozen=True)
class AnalysisConfig:
    event: str = "GW190521"
    gps: float = 1242442967.4
    ifos: tuple[str, ...] = ("H1", "L1", "V1")
    sample_rate: int = 2048
    event_duration: int = 32
    psd_duration: int = 512
    psd_gap: int = 32
    psd_fftlength: int = 8
    psd_overlap: int = 4
    f_low: float = 20.0
    f_high: float = 256.0
    on_source_half_width: float = 0.15
    redshift: float = 0.56
    mass_1_source: float = 98.4
    mass_2_source: float = 57.2
    distance_mpc: float = 3310.0
    # Illustrative fixed geometry, not an inferred sky position.
    right_ascension: float = 3.6
    declination: float = -0.7
    polarisation: float = 1.5
    inclination: float = 1.0
    burst_frequency: float = 60.0
    burst_q: float = 8.0
    burst_ellipticity: float = 1.0
    data_dir: Path = field(default_factory=lambda: Path("data"))
    figures_dir: Path = field(default_factory=lambda: Path("figures"))
    results_dir: Path = field(default_factory=lambda: Path("results"))

    @property
    def detector_mass_1(self):
        return self.mass_1_source * (1 + self.redshift)
    @property
    def detector_mass_2(self):
        return self.mass_2_source * (1 + self.redshift)
    @property
    def delta_t(self):
        return 1.0 / self.sample_rate
    @property
    def event_start(self):
        return int(self.gps) - self.event_duration // 2
    @property
    def psd_end(self):
        return int(self.gps) - self.psd_gap
    @property
    def psd_start(self):
        return self.psd_end - self.psd_duration
