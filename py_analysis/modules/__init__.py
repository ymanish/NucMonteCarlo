import sys
import os

sys.path.insert(0, os.path.expanduser("~/pol/Projects/Codebase/Spermatogensis/backend/NucFreeEnergy"))

NUC_STATE_PATH = os.path.expanduser("~/pol/Projects/Codebase/Spermatogensis/backend/NucFreeEnergy/methods/State/Nucleosome.state")
K_POSRESC_PATH = os.path.expanduser('~/pol/Projects/Codebase/Spermatogensis/backend/NucFreeEnergy/MDParams/nuc_K_posresc.npy')

__all__ = ['NUC_STATE_PATH', 'K_POSRESC_PATH']  

#### __all__ is a list of strings that explicitly declares which names (variables, functions, classes) should be considered "public" 
# when a user imports the module with from module import *.
#  If __all__ is defined, from module import * will only import the names listed in __all__.