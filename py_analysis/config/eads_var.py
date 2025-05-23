from typing import Optional
from dataclasses import dataclass, field

PHOSPHATE_SITES = [
    2, 6, 14, 17, 24, 29, 34, 38,
    45, 49, 55, 59, 65, 69, 76,
    80, 86, 90, 96, 100, 107, 111,
    116, 121, 128, 131, 139, 143,
]


@dataclass
class Parameters:
    NUCMETHOD: str
    FREEDNA_METHOD: Optional[str] = None 
    KRESCFACTOR: float= 1.0
    FLIP: bool = True
    HARD_CONS: bool = False

    TAU0: float = field(init=False) 

    def __post_init__(self):
        """Calculate TAU0 after instance initialization"""
        if self.NUCMETHOD == "hybrid":
            self.TAU0 = 1e-6
        elif self.NUCMETHOD == "crystal":
            self.TAU0 = 1e-6
        else:
            self.TAU0 = 1e-6

