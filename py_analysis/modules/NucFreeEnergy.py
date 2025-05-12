import os

# #Only import env_settings if the flag is set (default is to import)
if os.environ.get("IMPORT_ENV_SETTINGS", "1") == "1":
    from py_analysis.config.env_settings import *  # Triggers env_settings import

# from py_analysis.config.env_settings import *  # Triggers env_settings import
import numpy as np # Inherits the thread limits
from typing import Optional

from py_analysis.config.custom_types import FreeEnergyResult
from py_analysis.config.seq_var import *
from py_analysis.modules import NUC_STATE_PATH, K_POSRESC_PATH


from methods import nucleosome_free_energy, read_nucleosome_triads, GenStiffness, calculate_midstep_triads
from binding_model import binding_model_free_energy
from methods.PolyCG.polycg import cgnaplus_bps_params


class NucleosomeBreath:
    def __init__(self,  
                nuc_method='crystal', 
                hang_dna_method='md', hang_stiff=False, cgnaplus=False):
        """ 1. If hang_stiff is True, the hang_dna_method is used to calculate the stiffness of the free DNA.
            2. cgnaplus is not used here properly so better keep it False.
            3. Kentries varible is not used in the code, it has been replaced by the K_POSRESC matrix.
        """
        
        print(f"Using nucleosome method: {nuc_method}")
        print(f"Using hang DNA method: {hang_dna_method}")
        print(f"Using hang stiffness: {hang_stiff}")


        self.genstiff_nuc = GenStiffness(method=nuc_method)
        if hang_stiff:
            self.genstiff_hang = GenStiffness(method=hang_dna_method)
        self.triadfn = NUC_STATE_PATH
        self.nuctriads = read_nucleosome_triads(self.triadfn)

        self.fn = K_POSRESC_PATH
        # self.Kmat = np.load(self.fn)
        self.cgnaplus = cgnaplus


    def calculate_free_energy_soft(self, seq601:str, left:int, right:int, 
                                   id:Optional[str]=None, subid:Optional[str]=None,
                                     kresc_factor:float = 1, style:str="b_index")-> FreeEnergyResult:
        
        if self.cgnaplus:
            gs,stiff = cgnaplus_bps_params(seq601,group_split=True)
        else:
            stiff, gs = self.genstiff_nuc.gen_params(seq601, use_group=True)
        
        
        if style == "b_index":
            ## this converts the bound index to how many open phosphate sites required to be there
            ## Example: left=1, right=11, means l_open=2, r_open=4 so you need 2 open phosphate sites on left and 4 on right
            l_open = 2*left
            r_open  = 28-(2*right)-2

        elif style == "ph_index":
            ## Same as the bi style but directly at the phosphate sites level
            ## Example: left=0, right=27, means l_open=0, r_open=0
            ## Example: left=1, right=26, means l_open=1, r_open=1
            l_open = left
            r_open  = 28-(right)-1
           
        elif style == "open_sites":
            ## Here, we directly provide the number of open phosphate sites from left and right
            ## Example: left=0, right=27, means l_open=0, r_open=27
            l_open = left
            r_open  = right

        else:
            raise ValueError("Invalid style. Use 'b_index' or 'open_sites' or 'ph_index'.")
        
        K_resc = np.load(self.fn)

        # F_dict = soft_free_energy(gs,stiff,l_open,r_open, self.nuctriads, K_resc)
        nuc_mu0 = calculate_midstep_triads(triad_ids = self.select_phosphate_bind_sites(), 
                                           nucleosome_triads = self.nuctriads)
        

        F_dict = binding_model_free_energy(
            gs,
            stiff,    
            nuc_mu0,
            K_resc*kresc_factor,
            left_open=l_open,
            right_open=r_open,
            use_correction=True,
        )


        # Kentries = [1,1,1,10,10,10]
        # diags = np.concatenate([Kentries]*len(nuc_mu0))
        # K = np.diag(diags)

        # F_dict = binding_model_free_energy(
        #     gs,
        #     stiff,    
        #     nuc_mu0,
        #     K,
        #     left_open=l_open,
        #     right_open=r_open,
        #     use_correction=True,
        # )



        F601 = F_dict['F']
        F_entrop = F_dict['F_entropy']
        F_entalap = F_dict['F_enthalpy']
        F_free = F_dict['F_freedna']

        return FreeEnergyResult(F601, F_entrop, F_entalap, F_free, id, subid)
    

    def calculate_free_energy_hard(self, seq147:str, left:int, right:int, id:Optional[str]=None, subid:Optional[str]=None)-> FreeEnergyResult:
        stiff, gs = self.genstiff_nuc.gen_params(seq147, use_group=True)



        mid_array = self.select_phosphate_bind_sites(left, right)
        F_dict = nucleosome_free_energy(gs, stiff, mid_array, self.nuctriads, use_correction=True)


        F601 = F_dict['F']
        F_entrop = F_dict['F_entropy']
        F_entalap = F_dict['F_enthalpy']
        F_free = F_dict['F_freedna']



        # return F601, F_entrop, F_entalap, F_free, F_Diff
        return FreeEnergyResult(F601, F_entrop, F_entalap, F_free, id, subid)


        
    def select_phosphate_bind_sites(self, left=0, right=13):

        phosphate_bind_sites = [2, 6, 14, 17, 24, 29, 34, 38, 
                                    45, 49, 55, 59, 65, 69, 76, 
                                    80, 86, 90, 96, 100, 107, 111, 
                                    116, 121, 128, 131, 139, 143]
        
        return phosphate_bind_sites[left*2:(right*2)+2]
    

if __name__ == "__main__":
    
    # Example usage
    nuc_breath = NucleosomeBreath(nuc_method='crystal', hang_dna_method='md', hang_stiff=False, cgnaplus=False)
    result = nuc_breath.calculate_free_energy_soft(seq601="ACGTACGTACGTACGTACGTACGTACGT", left=1, right=2)
    print(result)