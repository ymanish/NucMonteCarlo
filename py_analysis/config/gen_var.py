
##### TETRAMER PARAMS #####
TETRAMER_PARAMS = {
     'NUC_PARAM_TYPE':"md", # "md" or "hybrid" or "crystal"
     'HANG_PARAM_TYPE':"md", # "md" or "hybrid" or "crystal"
     'TETRAMER_LENGTH':57,
     'NUC_LENGTH':147,
     'PAD_CHAR':"A",
     'BIND_POINTS':[4, 9], ## tetramer positions
     "HARD_CONS":False 
}    

# KENTRIES= np.array([1,1,1,10,10,10])*0.01
# KRESCFACTOR = 1.0

##### TI PARAMETERS #####
TI_PARAMS = {
     'SEQ_LENGTH': 147,
     'SLIDE_INTERVAL': 1,
     'WHOLE_SEQ': False,
     'MU_RANGE': [0, 1],
     'MU_VALUES': 20,
     'HARD_CONS': True
}
