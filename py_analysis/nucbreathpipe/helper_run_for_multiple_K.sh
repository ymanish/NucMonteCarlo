#### Helper script to run the nucleosome breath simulation for multiple K values in parallel in the local system.
# --- Conda Setup ---
CONDA_ENV="nucleosome"  
CONDA_SCRIPT_PATH="$HOME/anaconda3/etc/profile.d/conda.sh" 
PYTHON_PATH="/home/pol_schiessel/maya620d/anaconda3/envs/nucleosome/bin/python"
LOG_DIR="$HOME/pol/Projects/Codebase/NucMonteCarlo/logs"
# Source conda (required for activation in scripts)
if [ -f "$CONDA_SCRIPT_PATH" ]; then
    source "$CONDA_SCRIPT_PATH"
else
    echo "Error: conda.sh not found at $CONDA_SCRIPT_PATH"
    exit 1
fi



K_VALUES=(0.1 0.5 1.0 2.0 5.0 10.0) ## These are scaling factors for the histone core softness. Large values mean a hard core, small values mean a softer core.


# Run 2 jobs in parallel (each uses 5 cores -> 2×5=10 < 11)
parallel -j 2 --progress --bar --joblog $LOG_DIR/joblog.txt "conda run -n $CONDA_ENV && $PYTHON_PATH run_seqbatch_nucbreath.py {}" ::: "${K_VALUES[@]}"
# parallel -j 1 --progress --bar --joblog joblog.txt "conda run -n $CONDA_ENV && $PYTHON_PATH ../tetramer_fe.py {}" ::: "${K_VALUES[@]}"







