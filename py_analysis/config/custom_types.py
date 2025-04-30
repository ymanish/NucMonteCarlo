
from typing import NamedTuple, List, Optional, Tuple, Dict

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
