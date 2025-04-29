// Header file

#include "dnaMC2.h"
#include <time.h>

// Includes
#include <signal.h>
#include <unistd.h>
#include <stdbool.h>

// MC Sweep size
const int ENSSTEP = 147;

gsl_rng * r = NULL;

// START OF SIGNAL HANDLING CODE
//
// consider as black box
int sigint_received = 0;
int sigterm_received = 0;
int sigquit_received = 0;

void handle_sigint(int sig)  { sigint_received = 1;  }
void handle_sigterm(int sig) { sigterm_received = 1; }
void handle_sigquit(int sig) { sigquit_received = 1; }

static void setup_signal_handler(int sig, void (*handler)(  )) {
    #if _POSIX_VERSION > 198800L
    struct sigaction action;
    
    action.sa_handler = handler;
    sigemptyset(&(action.sa_mask));
    sigaddset(&(action.sa_mask), sig);
    action.sa_flags = 0;
    sigaction(sig, &action, 0);
    #else
    signal(sig, handler);
    #endif
}

// static int signal_was_caught(void)
// {
//     if (sigint_received) fprintf(stderr, "SIGINT received!\n");
//     if (sigterm_received) fprintf(stderr, "SIGTERM received!\n");
//     if (sigquit_received) fprintf(stderr, "SIGQUIT received!\n");
//     return (sigint_received || sigterm_received || sigquit_received);
// }
// END OF SIGNAL HANDLING CODE



// Globally accessible sequence arrays
Base *Sequence;
Base *PullSequence;
Base *Sequence2;
Base *PullSequence2;

// Array that keeps track of the right-hand base pairs of the 28 base pair steps that have
// fixed midplanes in the RBP nucleosome model.
// int boundbp[29] = {5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000};

int boundbp[29] = {3, 7, 15, 18, 25, 30, 35, 39, 46, 50, 56, 60, 66, 70, 77, 81, 87, 91, 97, 101, 108, 112, 117, 122, 129, 132, 140, 144, 5000};



// Initialization function. Currently only sets up the 
// random number generator.
void init(int seed) {
    
    gsl_rng_env_setup();
    r = gsl_rng_alloc (gsl_rng_mt19937);
    gsl_rng_set (r, seed);
    
}


// Signal handling function
void signal_handler(int signal_num) {
    printf("Interrupt signal is (%d).\n", signal_num);
    exit(signal_num);  
}




int main(int argc, char *argv[]){
    
    signal(SIGINT, signal_handler);
    

    // DNAState arrays
    DNAState state[N];
    DNAState start[N];
    DNAState trial[N];
    DNAState midstate[N-1];
    DNAState midplanes[28];
    
    double BPSEnergies[N-1], BPSEnergiesTrial[N-1];
    double BPSEnergies2[N-1], BPSEnergiesTrial2[N-1];

    
    int i, n, n1, l, id_S, id_E, sweep_size, LED_n1, RED_n1;
    int bb, bplo, bphi;
    int seed = 10089; 
    int startpos, nstep;
    int base_sel_left, base_sel_right;

    double j;
    double beta = 1;

    double stepsize_t = 1.0/(sqrt(beta)*10.0); 
    double stepsize_a = PI/(sqrt(beta)*300.0);

    double e_A, e_B, d, w, dE_A, dE_B, dH, s, H_e;
    
    // MC loop counter and length
    long int sweeps_c, midbp_tries, sbp_tries, tries_crank, tries_cluster, midbp_accept, sbp_accept, accept_crank, accept_cluster, equi_sweeps;
    long int stop = 1e6;

    double diff, mean_energy, variance, std_err;


    // Auxiliary variables for reading in sequences

    int numCharacters = 0, cnt = 0;
    double p;
     
   
    static double sd_energy_sum, sd_energy_sum_squared, energy_diff;
    static long int sd_energy_samples, sampled_data;
    static bool equ_flag, break_flag;





    ////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////

    
    if (argc > 1) {
        startpos = atoi(argv[1]);
        nstep = atoi(argv[2]);

        if (argc > 3) {
            base_sel_left = atoi(argv[3]);
            base_sel_right = atoi(argv[4]);
        }
        if (argc > 4) {
            p = atof(argv[5]);
        }
    }
    else { // Default is to start at the beginning and take only the first position
        startpos = 0;
        nstep = 1;
        base_sel_left = 0;
        base_sel_right = 13;
        p = 0.5;

    }
    
    compute_bound_sites(base_sel_left, base_sel_right);

    // Debug: Print the modified boundbp array
    printf("Modified boundbp array:\n");
    for (int i = 0; i < 29; i++) {
        printf("%d ", boundbp[i]);
    }
    printf("\n");

   
    Sequence += startpos;
    Sequence2  += startpos;

    // Initialization routines
    init(seed);

    
    fillK();


    // Read in sequences from SequenceA and SequenceB files
    const char *filenames[] = {"./State/SequenceA.seq", "./State/SequenceB.seq"};
    Base **sequences[] = {&Sequence, &Sequence2};
    Base **pullSequences[] = {&PullSequence, &PullSequence2};

    for (int fileIndex = 0; fileIndex < 2; fileIndex++) {
        numCharacters = 0;
        cnt = 0;

        // Read in sequence length
        fprintf(stderr, "Reading sequence length from %s.\n", filenames[fileIndex]);
        FILE* seqfile = fopen(filenames[fileIndex], "r");
        if (seqfile == NULL) {
            fprintf(stderr, "Error: Could not open file %s for reading.\n", filenames[fileIndex]);
            exit(1);
        }
        char nextChar = getc(seqfile);
        while (nextChar != EOF) {
            if (isspace(nextChar)) {
                nextChar = getc(seqfile);
            } else {
                numCharacters++;
                nextChar = getc(seqfile);
            }
        }
        fclose(seqfile);

        // Allocate memory to hold the sequence
        fprintf(stderr, "Allocating space for sequence from %s.\n", filenames[fileIndex]);
        *pullSequences[fileIndex] = (Base*) malloc(numCharacters * sizeof(int));
        *sequences[fileIndex] = *pullSequences[fileIndex];

        // Read in the sequence again and store in memory
        fprintf(stderr, "Reading sequence from %s.\n", filenames[fileIndex]);
        seqfile = fopen(filenames[fileIndex], "r");
        if (seqfile == NULL) {
            fprintf(stderr, "Error: Could not open file %s for reading.\n", filenames[fileIndex]);
            exit(1);
        }
        nextChar = getc(seqfile);
        while (nextChar != EOF) {
            if (isspace(nextChar)) {
                nextChar = getc(seqfile);
            } else {
                (*sequences[fileIndex])[cnt] = CharToBase(nextChar);
                cnt++;
                nextChar = getc(seqfile);
            }
        }
        fclose(seqfile);
        fprintf(stderr, "Sequence read from %s.\n", filenames[fileIndex]);
    }
    

    fprintf(stderr, "Reading state.\n");
    FILE* statefile = fopen("State/Nucleosome.state", "r");



    PrintSequence(stdout, Sequence);
    PrintSequence(stdout, Sequence2);


    ReadDNAState(statefile, state);



    fclose(statefile);
    fprintf(stderr, "State read.\n");
    


    // Calculate midplanes for the state just read in
    midplanestate(state, midstate);

    
    // Store the 28 relevant midframes in a new array
    for ( i = 0; i < N; i++ ) {
        if ( BaseToBond(i) >= 0 ) {
            //fprintf(stderr, "%d\n", i);
            midplanes[BaseToBond(i)] = midstate[i-1];
        }
    }
    





    ////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////



    fprintf(stdout, "Starting MC.\n");

    // Initialize variables
    fprintf(stdout, "Calculating for p = %.3f\n", p);


    char energy_filename[200];
    snprintf(energy_filename, sizeof(energy_filename), "./TI_output/Hamiltonian_Energies_p_%.3f.dat", p);
    FILE *energy_file = fopen(energy_filename, "w");
    if (energy_file == NULL) {
        fprintf(stderr, "Error: Could not open energies.dat for writing.\n");
        exit(1);
    }
    
    for ( n = 0; n < N; n++ ) {
        // For safety, orthogonalize the orientations before starting
        OrthogonalizeR(state[n].R);
    }

    CopyDNAState(state, trial, 0, N-1);
    
    PrintBPState_terminal(state[0]);  // Print the contents of myState
    PrintBPState_terminal(state[37]); 
    
    // Start energy tracking
    e_A = Energy(state, Sequence);
    e_B = Energy(state, Sequence2);
    H_e = (1-p)*e_A + p*e_B;


    fprintf(stdout, "Initial Energy: %f\t%f\t%f\n", H_e, e_A, e_B);

    for ( n = 0; n < N-1; n++ ) {
       
        BPSEnergies[n] = BPSEnergy(state[n], state[n+1], Sequence[n], Sequence[n+1]);
        BPSEnergiesTrial[n] = BPSEnergies[n];


        BPSEnergies2[n] = BPSEnergy(state[n], state[n+1], Sequence2[n], Sequence2[n+1]);
        BPSEnergiesTrial2[n] = BPSEnergies2[n];

    }


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    sweep_size = 1000; //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    equi_sweeps= 10000;
    
    sbp_tries = 0; //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    midbp_tries = 0;

    // tries = 0;
    sbp_accept = 0; //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    midbp_accept = 0;

    tries_crank = 0;
    tries_cluster = 0; //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    accept_crank = 0;
    accept_cluster = 0; //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    energy_diff=0.0; //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    sampled_data=0;
    equ_flag=false;
    break_flag=false;

    dE_A = 0.0; //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    dE_B = 0.0;
    dH = 0.0;
    sweeps_c=0;
    std_err = 100.0;

    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    fprintf(stderr, "Starting Main LOOP for MC.........\n");

    while (std_err > 0.01) {    
        
        sweeps_c += 1;
        for (int move = 0; move < sweep_size; move++) {

            
            for (int bp_move = 0; bp_move < 147; bp_move++) {
                
                n1 = gsl_rng_uniform_int(r, N);
                perform_sbp_move(n1, stepsize_a, stepsize_t, beta,
                        state,
                        trial,
                        BPSEnergies,
                        BPSEnergiesTrial,
                        BPSEnergies2,
                        BPSEnergiesTrial2,
                        midplanes,
                        p,
                        &e_A,
                        &e_B,
                        &H_e,
                        &sbp_tries,
                        &midbp_tries,
                        &sbp_accept,
                        &midbp_accept); 
            }


            LED_n1 = gsl_rng_uniform_int(r, boundbp[0]);

            perform_sbp_move(LED_n1, stepsize_a, stepsize_t, beta,
                            state,
                            trial,
                            BPSEnergies,
                            BPSEnergiesTrial,
                            BPSEnergies2,
                            BPSEnergiesTrial2,
                            midplanes,
                            p,
                            &e_A,
                            &e_B,
                            &H_e,
                            &sbp_tries,
                            &midbp_tries,
                            &sbp_accept,
                            &midbp_accept); 

            RED_n1 = gsl_rng_uniform_int(r, N-boundbp[28]) + boundbp[28];
            perform_sbp_move(RED_n1, stepsize_a, stepsize_t, beta,
                            state,
                            trial,
                            BPSEnergies,
                            BPSEnergiesTrial,
                            BPSEnergies2,
                            BPSEnergiesTrial2,
                            midplanes,
                            p,
                            &e_A,
                            &e_B,
                            &H_e,
                            &sbp_tries,
                            &midbp_tries,
                            &sbp_accept,
                            &midbp_accept); 


            ///////////////////////////////////////////////////////////////////////
            ///////////////////////////////////////////////////////////////////////
            ////////////////////////CRANKSHAFT MOVE AND CLUSTER TRANS /////////////
            ///////////////////////////////////////////////////////////////////////
        


            performCrankshaftMoves(stepsize_a, beta,
                                    state, trial,
                                    BPSEnergies,  BPSEnergiesTrial,
                                    BPSEnergies2, BPSEnergiesTrial2,
                                    p, &e_A, &e_B, &H_e,
                                    &tries_crank, &accept_crank);



            perform_cluster_trans_moves(stepsize_t, beta, p,
                                    state, trial,
                                    BPSEnergies, BPSEnergiesTrial,
                                    BPSEnergies2, BPSEnergiesTrial2,
                                    &e_A, &e_B, &H_e,
                                    &tries_cluster, &accept_cluster);

        
        
        
        }


        // ONE SWEEP DONE

        // fprintf(stdout, "Sweep Number: %ld\t H: %e\t E_A: %e\t E_B: %e\n", sweeps_c, H_e, e_A, e_B);


        if (equ_flag==true && sweeps_c % 1 == 0){

            diff = e_B-e_A;
            energy_diff += diff;
            // sampled_data ++;

            sd_energy_sum += H_e;
            sd_energy_sum_squared += H_e * H_e;
            sd_energy_samples++;

            mean_energy = sd_energy_sum / sd_energy_samples;
            variance = (sd_energy_sum_squared / sd_energy_samples) - (mean_energy * mean_energy);
            std_err = sqrt(variance) / sqrt(sd_energy_samples);

            // fprintf(energy_file, "%.3f,%.6f\n", p, mean_energy);
            // fprintf(energy_file, "%ld,%.6f\n", sweeps_c, H_e);
            
        }



        else if (sweeps_c > equi_sweeps && equ_flag==false) {
            equ_flag = true;
            // fprintf(stdout, "Equilibration achieved. Mean: %.6f, Standard Error: %.6f\n", mean_energy, std_err);
            fprintf(stdout, "Equilibration achieved\n");
            fprintf(stdout, "Sweep Number: %ld\t H: %e\t E_A: %e\t E_B: %e\n", sweeps_c, H_e, e_A, e_B);
            fprintf(stdout, "Sampling Data...\n");
            sd_energy_sum = H_e;
            sd_energy_sum_squared = H_e * H_e;
            sd_energy_samples = 1;
            }

        // else {
        //     fprintf(stdout, "Equilibrating...\n");
        //     }

        
        if (sweeps_c % 1000 == 0){
            time_t now = time(NULL);
            char *timestamp = ctime(&now);
            timestamp[strlen(timestamp) - 1] = '\0'; // Remove the newline character
            // double rate = (double)accept / (double)tries;
            double sbp_rate = (double)sbp_accept / (double)sbp_tries;
            double mdp_rate = (double)midbp_accept / (double)midbp_tries;
            double crank_rate = (double)accept_crank / (double)tries_crank;
            double cluster_rate = (double)accept_cluster / (double)tries_cluster;

            fprintf(stdout, "[%s] Sweep Number:%ld, Mean: %.6f, Variance: %.6f, Standard Error: %.6f\n", timestamp, sweeps_c, mean_energy, variance, std_err);
            fprintf(stdout, "Acceptance Rates.................\n");
            fprintf(stdout, "SBP_Rate: %.2f, MIDBase_Rate:%.2f, CrankShaft_rate: %.2f, ClusterTrans_rate: %.2f\n\n", sbp_rate, mdp_rate, crank_rate, cluster_rate);

            // fprintf(stderr, "-------------------------------------\n");
            }

        
        // If any signal was caught, break out of the loop
        if (signal_was_caught(  )) break;

        if (sweeps_c > 1e6) {
            fprintf(energy_file, "%ld,%.6f\n", sweeps_c, H_e);
            fprintf(stdout, "End of simulation: %ld\n", sweeps_c);
            return 0;
        }

        
    
    }

    mean_energy = energy_diff / sd_energy_samples;
    fprintf(energy_file, "%.3f,%.6f\n", p, mean_energy);
    fprintf(stdout, "Mean energy difference: %.6f\n", mean_energy);
    fprintf(stdout, "End of simulation: %.3f\n", p);
    
    // combineXYZFiles();
    // Close the file
    fclose(energy_file);
    
    // Clean up the random number generator.
    gsl_rng_free (r);
    
    // Cleanly exit, returning exit status 0 (success).
    return 0;
}

