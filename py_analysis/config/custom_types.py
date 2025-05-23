
from typing import NamedTuple, Optional
import pandas as pd
from dataclasses import dataclass, field
from enum import Enum, auto
import numpy as np
import math

class FreeEnergyResult(NamedTuple):
    """Container for energy calculation results."""
    F: float
    F_entropy: float
    F_enthalpy: float
    F_freedna: float
    id: Optional[str]
    subid: Optional[str]

class ProcessedSequence(NamedTuple):
    """Container for a processed sequence and its metadata."""
    id: str
    subid: str
    sequence: str
    start_site: int
    end_site: int

class NuclBreathingResult(NamedTuple):
    """Container for nucleosome breathing results."""
    id: str
    subid: str
    sequence: str
    leftbind_indx: int
    rightbind_indx: int
    F_vals: FreeEnergyResult
    Adsorp_F: float = 0.0

##################################################
##################################################
########## DWELL TIME: GAUSSIAN FITTING ##########
##################################################
##################################################
##################################################

### slots=True:
### Instances use __slots__ instead of __dict__.
### You can modify declared fields, but cannot add new ones.
### Memory per instance is smaller; attribute lookup is slightly faster.
@dataclass(slots=True)
class FitResult:
    """Holds outputs from *fit_multi_gaussian*."""

    centres: np.ndarray        # bp positions used as initial centres
    pars: np.ndarray           # raw parameters from curve_fit
    pcov: np.ndarray           # covariance matrix

    @property
    def amplitudes(self) -> np.ndarray:  # A
        return self.pars[0::3]

    @property
    def mus(self) -> np.ndarray:  # fitted centres (mu)
        return self.pars[1::3]

    @property
    def sigmas(self) -> np.ndarray:  # sigma
        return self.pars[2::3]

    @property
    def areas(self) -> np.ndarray:
        """Analytic integral under each Gaussian peak."""
        return self.amplitudes * self.sigmas * math.sqrt(2.0 * math.pi)

    # convenience for tabular view
    @property
    def results(self) -> pd.DataFrame:
        return pd.DataFrame({
            "dwell_time": self.areas,
            "peak_pos":   np.round(self.mus, 0),
            "sigma":      self.sigmas,
            "A":          self.amplitudes,
        })


@dataclass(slots=True)
class DwellTimeResult:
    """Container for final results for the dwell times after gaussian fitting."""
    phosphate_sites: np.ndarray
    left_dwell_times: np.ndarray
    right_dwell_times: np.ndarray
    
    def to_dataframe(self) -> pd.DataFrame:
        return pd.DataFrame({
            "phosphate_site": self.phosphate_sites,
            "left_dwell_time": self.left_dwell_times,
            "right_dwell_time": self.right_dwell_times
        })


##################################################
##################################################
########## ADSORPTION ENERGY CALCULATION #########
##################################################
##################################################
##################################################

class Orientation(Enum):
    """Physical orientation of the DNA on the histone octamer."""

    ORIGINAL = auto()  # -73 -> +73 (as in the crystal structure)
    FLIPPED = auto()   # same octamer rotated 180 degrees, so +73 -> -73

    @property
    def is_flipped(self) -> bool:
        return self is Orientation.FLIPPED

#### frozen=True: Instances remain dict-based but become immutable: you cannot change or add fields once constructed.
@dataclass(frozen=True)
class EadsResult:
    """Computed adsorption energies for *one* orientation."""

    orientation: Orientation

    Eads_Base_left: pd.Series
    Eads_Base_right: pd.Series

    Eads_dG_left: pd.Series
    Eads_dG_right: pd.Series

    Eads_dH_left: pd.Series
    Eads_dH_right: pd.Series

    # raw pieces (handy for debugging / further analysis)
    tau_left: pd.Series = field(repr=False)
    tau_right: pd.Series = field(repr=False)
    dG_left: pd.Series = field(repr=False)
    dG_right: pd.Series = field(repr=False)
    dH_left: pd.Series = field(repr=False)
    dH_right: pd.Series = field(repr=False)

