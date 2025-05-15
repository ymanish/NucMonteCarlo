from dataclasses import dataclass, field
from pathlib import Path
import datetime
from typing import ClassVar, Dict
######LOCAL DIR#######

SRC_DIR = Path(__file__).parent.parent.parent
OUT_DIR = Path(__file__).parent.parent.parent

######CLUSTER DIR#######

# SRC_DIR = Path('/home/pol_schiessel/maya620d/NucMonteCarlo')
# OUT_DIR = Path('/group/pol_schiessel/Manish/NucMonteCarlo')


print(f"Source Directory: {SRC_DIR}")
print(f"Output Directory: {OUT_DIR}")


BACKEND_DIR = SRC_DIR / 'backend'
PARAM_DIR = SRC_DIR / 'parameters'
TEST_DIR = SRC_DIR / 'tests'
EXAMPLES_DIR = SRC_DIR / 'examples'


DATA_DIR = OUT_DIR / 'data'
RESULTS_DIR = OUT_DIR / 'output'




##### FILE NAMES #####


@dataclass(frozen=True)
class ddGDataSaveParams:
    nuc_type: str
    hangdna_type: str
    tetramer_length: int
    pad_char: str
    constraint: str
    kresc_factor: float

     # Class-level configuration
    ABBREVIATIONS: ClassVar[Dict[str, str]] = field(
        default={
            'constraint': 'c',
            'nuc': 'n',
            'freedna': 'fd',
            'tet': 't',
            'pad': 'p',
            'scatter': 'sct',
            '_': '-',  # Separator, 
            'kfactor': 'kf',
        },
        repr=False
    )


    
    def ddG_filename(self, long_name: bool = False) -> str:

        if long_name:
            parts = [f"scatter_nuc_{self.nuc_type}_"
                f"freedna_{self.hangdna_type}_"
                f"constraint_{self.constraint}_"
                f"tet_{self.tetramer_length}_"
                f"pad_{self.pad_char}_"
                f"kfactor_{self.kresc_factor}"]
                
            
        else:

            parts = [
                f"{self.ABBREVIATIONS['constraint']}{self.constraint}",
                f"{self.ABBREVIATIONS['nuc']}{self.nuc_type}",
                f"{self.ABBREVIATIONS['freedna']}{self.hangdna_type}",
                f"{self.ABBREVIATIONS['tet']}{self.tetramer_length}",
                f"{self.ABBREVIATIONS['pad']}{self.pad_char}",
                f"{self.ABBREVIATIONS['kfactor']}{self.kresc_factor}",
                ]
        return "_".join(parts)

  
