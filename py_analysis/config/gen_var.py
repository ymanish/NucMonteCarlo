from dataclasses import dataclass, field
from pathlib import Path
import datetime
from typing import ClassVar, Dict

NUC_PARAM_TYPE = "hybrid" # "md" or "hybrid" or "crystal"
HANG_PARAM_TYPE = "md" # "md" or "hybrid" or "crystal"
TETRAMER_LENGTH = 57
NUC_LENGTH = 147
PAD_CHAR = "A"
BIND_POINTS = [4, 9] ## tetramer positions

# KENTRIES= np.array([1,1,1,10,10,10])*0.01

select_min_F = False
HARD_CONS = False


##### TI PARAMETERS #####

TI_PARAMS = {
     'SEQ_LENGTH': 147,
     'SLIDE_INTERVAL': 1,
     'WHOLE_SEQ': False,
     'MU_RANGE': [0, 1],
     'MU_VALUES': 20,
}

######LOCAL DIR#######

SRC_DIR = Path(__file__).parent.parent.parent
OUT_DIR = Path(__file__).parent.parent.parent

######CLUSTER DIR#######

# SRC_DIR = Path('/home/pol_schiessel/maya620d/')
# OUT_DIR = Path('/group/pol_schiessel/Manish/')


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

     # Class-level configuration
    ABBREVIATIONS: ClassVar[Dict[str, str]] = field(
        default={
            'constraint': 'c',
            'nuc': 'n',
            'freedna': 'fd',
            'tet': 't',
            'pad': 'p',
            'scatter': 'sct',
            '_': '-'  # Separator
        },
        repr=False
    )


    
    def ddG_filename(self, long_name: bool = False) -> str:

        if long_name:
            parts = [f"scatter_nuc_{self.nuc_type}_"
                f"freedna_{self.hangdna_type}_"
                f"constraint_{self.constraint}_"
                f"tet_{self.tetramer_length}_"
                f"pad_{self.pad_char}"]
                
            
        else:

            parts = [
                f"{self.ABBREVIATIONS['constraint']}{self.constraint}",
                f"{self.ABBREVIATIONS['nuc']}{self.nuc_type}",
                f"{self.ABBREVIATIONS['freedna']}{self.hangdna_type}",
                f"{self.ABBREVIATIONS['tet']}{self.tetramer_length}",
                f"{self.ABBREVIATIONS['pad']}{self.pad_char}"
                ]
        return "_".join(parts)

  