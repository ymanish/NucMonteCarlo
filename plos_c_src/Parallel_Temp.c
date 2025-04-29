// Header file

#include "dnaMC2.h"
#include <sys/time.h>
#include <signal.h>
#include <unistd.h>
#include <stdbool.h>
#include <omp.h>

// MC Sweep size
const int ENSSTEP = 147;

gsl_rng * rg = NULL;
// gsl_rng * r2 = NULL;

// double betas[NUM_REPLICAS] = {1.0,0.9,0.8,0.7};

// double Mus[NUM_REPLICAS] = {0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0};
double Mus[NUM_REPLICAS] = {0.5};

// double Mus[NUM_REPLICAS] = {0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};
double energies[NUM_REPLICAS] = {0.0};

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

int signal_was_caught(void)
{
    if (sigint_received) fprintf(stderr, "SIGINT received!\n");
    if (sigterm_received) fprintf(stderr, "SIGTERM received!\n");
    if (sigquit_received) fprintf(stderr, "SIGQUIT received!\n");
    return (sigint_received || sigterm_received || sigquit_received);
}


// Globally accessible sequence arrays
Base *Sequence;
Base *PullSequence;
Base *Sequence2;
Base *PullSequence2;

int global_boundbp[29] = {3, 7, 15, 18, 25, 30, 35, 39, 46, 50, 56, 60, 66, 70, 77, 81, 87, 91, 97, 101, 108, 112, 117, 122, 129, 132, 140, 144, 5000};

Replica *replicas[NUM_REPLICAS];

// Initialization function. Currently only sets up the 
// random number generator.
void init(int seed) {
    
    gsl_rng_env_setup();
    rg = gsl_rng_alloc (gsl_rng_mt19937);
    gsl_rng_set (rg, seed);
    
}

// Signal handling function
void signal_handler(int signal_num) {
    printf("Interrupt signal is (%d).\n", signal_num);
    exit(signal_num);  
}


static void freeReplica(Replica *rep) {
    if (rep == NULL) return;

    // Free dynamically allocated arrays
    free(rep->state);
    free(rep->trial);
    free(rep->midstate);
    free(rep->midplanes);
    free(rep->BPSEnergies);
    free(rep->BPSEnergies2);
    free(rep->BPSEnergiesTrial);
    free(rep->BPSEnergiesTrial2);
    free(rep->boundbp);
    free(rep->rn);

    // Free the Replica object itself
    free(rep);
}


Replica* createReplica() {
    Replica *rep = malloc(sizeof(Replica));  // Allocate a Replica object on the heap.
    if (!rep) { perror("malloc"); exit(1); }
    
    // Allocate memory for each of the arrays inside the replica.
    rep->beta = 1.0;
    rep->state = malloc(N * sizeof(DNAState));
    rep->trial = malloc(N * sizeof(DNAState));
    rep->midstate = malloc((N - 1) * sizeof(DNAState));
    rep->midplanes = malloc(28 * sizeof(DNAState));
    rep->BPSEnergies = malloc((N - 1) * sizeof(double));
    rep->BPSEnergies2 = malloc((N - 1) * sizeof(double));
    rep->BPSEnergiesTrial = malloc((N - 1) * sizeof(double));
    rep->BPSEnergiesTrial2 = malloc((N - 1) * sizeof(double));
    rep->energy_A = 0.0;
    rep->energy_B = 0.0;
    rep->H_e = 0.0;
    rep->swap_tries = 0;
    rep->accepted_swaps = 0;
    rep->equilibrated = 0;

    rep->energy_stats.sample=0;
    rep->energy_stats.value_sum=0;
    rep->energy_stats.sample_mean=0;
    rep->energy_stats.sum_sq_diff=0;
    rep->energy_stats.sample_variance=0;
    rep->energy_stats.sample_std_dev=0;
    rep->energy_stats.standard_error=0;
    rep->mu = 0.0;
    rep->rn = NULL;

    rep->boundbp = malloc((29) * sizeof(int));
   

    return rep;
}


void TEMPcopyState(DNAState *dest, const DNAState *src, int L) {
    for (int i = 0; i < L; i++) {
        dest[i] = src[i];
    }
}

void TEMPcopyEnergies(double *dest, const double *src, int L) {
    for (int i = 0; i < L; i++) {
        dest[i] = src[i];
    }
}

void update_running_stats(RunningStats *stats, double new_value) {
    stats->sample++;
    stats->value_sum+= new_value;
    stats->sample_mean = stats->value_sum / stats->sample;
    stats->sum_sq_diff += new_value * new_value;


    stats->sample_variance = (stats->sum_sq_diff/stats->sample)-(stats->sample_mean*stats->sample_mean);
    stats->sample_std_dev = sqrt(stats->sample_variance);
    stats->standard_error = stats->sample_std_dev / sqrt(stats->sample);

}


int is_equilibrated(const Replica *replica, double rel_tol) {
    if (replica->energy_stats.sample < 1000) return 0;  // Minimum samples

    return (replica->energy_stats.sample_std_dev / fabs(replica->energy_stats.sample_mean) < rel_tol) ? 1:0;  // e.g., rel_tol = 0.01 (1%)
}

void printReplicas(Replica **replicas, int num_replicas) {
    for (int i = 0; i < num_replicas; i++) {
        printf("Replica %d Details:\n", i);
        printf("  Beta            : %f\n", replicas[i]->beta);
        printf("-----Mu-----      : %f\n", replicas[i]->mu);
        printf("  Boundbp         : ");
        for (int j = 0; j < 29; j++) {
            printf("%d ", replicas[i]->boundbp[j]);
        }
        printf("\n");
        printf("  Energy_A        : %f\n", replicas[i]->energy_A);
        printf("  Energy_B        : %f\n", replicas[i]->energy_B);
        printf("  H_e             : %f\n", replicas[i]->H_e);
        printf("  Accepted Swaps  : %ld\n", replicas[i]->accepted_swaps);
        printf("  Energy Stats:\n");
        printf("    Samples           : %ld\n", replicas[i]->energy_stats.sample);
        printf("    Sample Mean       : %f\n", replicas[i]->energy_stats.sample_mean);
        printf("    Standard Deviation: %f\n", replicas[i]->energy_stats.sample_std_dev);
        printf("    Standard Error    : %f\n", replicas[i]->energy_stats.standard_error);
        printf("\n");
    }

}



int main(int argc, char *argv[]){

    struct timeval tv_start, tv_end;
    gettimeofday(&tv_start, NULL);

    signal(SIGINT, signal_handler);

     // DNAState arrays
    DNAState global_state[N];
    DNAState global_trial[N];
    DNAState global_midstate[N-1];
    DNAState global_midplanes[28];
    
    double global_BPSEnergies[N-1], global_BPSEnergiesTrial[N-1];
    double global_BPSEnergies2[N-1], global_BPSEnergiesTrial2[N-1];
    
    int i, n, n1, l, id_S, id_E;
    int seed = 10089; 
    int startpos, nstep;
    int base_sel_left, base_sel_right;
    double beta;
    double e_A, e_B, H_e;
    long int TOT_SWEEP = 1e6;
    long int EQUI_SWEEP = 1e5;
    double rel_tol = 0.01;
    int numCharacters = 0, cnt = 0;

    long int tries = 0;

     
    long int sbp_tries = 0; 
    long int midbp_tries = 0;

    long int sbp_accept = 0; 
    long int midbp_accept = 0;

    long int tries_crank = 0;
    long int tries_cluster = 0; 

    long int accept_crank = 0;
    long int accept_cluster = 0; 




    if (argc > 1) {
        startpos = atoi(argv[1]);
        nstep = atoi(argv[2]);

        if (argc > 3) {
            base_sel_left = atoi(argv[3]);
            base_sel_right = atoi(argv[4]);
        }
        // if (argc > 4) {
        //     p = atof(argv[5]);
        // }
    }
    else { // Default is to start at the beginning and take only the first position
        startpos = 0;
        nstep = 1;
        base_sel_left = 0;
        base_sel_right = 13;

    }
    

    // // Debug: Print the modified boundbp array
    // printf("Modified boundbp array:\n");
    // for (int i = 0; i < 29; i++) {
    //     printf("%d ", boundbp[i]);
    // }
    // printf("\n");

   
    Sequence += startpos;
    Sequence2  += startpos;

    // Initialization routines
    init(seed);
    // init(seed+1);

    
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


    ReadDNAState(statefile, global_state);



    fclose(statefile);
    fprintf(stderr, "State read.\n");
 
    // Calculate midplanes for the state just read in
    midplanestate(global_state, global_midstate);

    
    // Store the 28 relevant midframes in a new array
    for ( int i = 0; i < N; i++ ) {
        if ( BaseToBond(i, global_boundbp) >= 0 ) {
            //fprintf(stderr, "%d\n", i);
            global_midplanes[BaseToBond(i, global_boundbp)] = global_midstate[i-1];
        }
    }
    


    for (int n = 0; n < N; n++ ) {
        // For safety, orthogonalize the orientations before starting
        OrthogonalizeR(global_state[n].R);
    }

    CopyDNAState(global_state, global_trial, 0, N-1);
    
    PrintBPState_terminal(global_state[0]);  // Print the contents of myState
    PrintBPState_terminal(global_state[37]); 
    
    // Start energy tracking
    e_A = Energy(global_state, Sequence);
    e_B = Energy(global_state, Sequence2);
    // H_e = (1-p)*e_A + p*e_B;


    fprintf(stdout, "Initial Energy for E_A and E_B: %f\t%f\n", e_A, e_B);

    for (int n = 0; n < N-1; n++ ) {
       
        global_BPSEnergies[n] = BPSEnergy(global_state[n], global_state[n+1], Sequence[n], Sequence[n+1]);
        global_BPSEnergiesTrial[n] = global_BPSEnergies[n];


        global_BPSEnergies2[n] = BPSEnergy(global_state[n], global_state[n+1], Sequence2[n], Sequence2[n+1]);
        global_BPSEnergiesTrial2[n] = global_BPSEnergies2[n];

    }



    for (int i = 0; i < NUM_REPLICAS; i++) {


        replicas[i] = createReplica();

        // replicas[i]->beta = betas[i];
        // replicas[i]->Sequence = Sequence;
        // replicas[i]->Sequence2 = Sequence2;
        replicas[i]->mu = Mus[i];

        TEMPcopyState(replicas[i]->state, global_state, N);
        TEMPcopyState(replicas[i]->trial, global_trial, N);
        TEMPcopyState(replicas[i]->midstate, global_midstate, N-1);
        TEMPcopyState(replicas[i]->midplanes, global_midplanes, 28);

        TEMPcopyEnergies(replicas[i]->BPSEnergies, global_BPSEnergies, N-1);
        TEMPcopyEnergies(replicas[i]->BPSEnergies2, global_BPSEnergies2, N-1);
        TEMPcopyEnergies(replicas[i]->BPSEnergiesTrial, global_BPSEnergiesTrial, N-1);
        TEMPcopyEnergies(replicas[i]->BPSEnergiesTrial2, global_BPSEnergiesTrial2, N-1);

        replicas[i]->energy_A = e_A;
        replicas[i]->energy_B = e_B;
        replicas[i]->H_e = (1-replicas[i]->mu)*e_A + replicas[i]->mu*e_B;

        // replicas[i]->accepted_swaps = 0;
        // replicas[i]->equilibriated = 0;
        // replicas[i]->energy_stats.sample=0;
        // replicas[i]->energy_stats.sample_mean=0;
        // replicas[i]->energy_stats.sum_sq_diff=0;


         // Allocate a unique RNG for each replica.
        replicas[i]->rn = gsl_rng_alloc(gsl_rng_mt19937);
        gsl_rng_set(replicas[i]->rn, seed + i);

        compute_bound_sites(base_sel_left, base_sel_right, replicas[i]->boundbp);


    }


    // fprintf(stdout, "Calculating for p = %.3f\n", p);



    // printReplicas(replicas, NUM_REPLICAS); 
    // exit(0);



//  #pragma omp parallel  // Creates a parallel region once
//     {
//         // #pragma omp single
//         // {
//         //     int num_threads = omp_get_num_threads();
//         //     printf("Number of threads activated for omp parallel: %d\n", num_threads);
//         // }
//         // printf("Hellooooooo\n");

//         for (int S=0; S<10; S++) {
            
//             // Code here runs only on one thread (typically thread 0)
//             #pragma omp single
//             {printf("[Sweep %d] Starting replica updates on thread %d\n", S, omp_get_thread_num());}
        

//             #pragma omp for schedule(dynamic)
//             for (int i = 0; i < 10; i++){
//                 printf("[Sweep %d] [Replica %d] running on thread %d\n", S, i, omp_get_thread_num());

//                 }

//             #pragma omp single
//             {printf("Done...............\n");}

//         }
        
//     }
//     exit(0);

    char energy_filename[200];
    snprintf(energy_filename, sizeof(energy_filename), "./TI_output/Auto_Correlation_Replica1_S100_sweeps10.dat");
    FILE *energy_file = fopen(energy_filename, "w");
    if (energy_file == NULL) {
        fprintf(stderr, "Error: Could not open energies.dat for writing.\n");
        exit(1);
    }

    char energy_filename1[200];
    snprintf(energy_filename1, sizeof(energy_filename1), "./TI_output/Auto_Correlation_Replica1_S0_sweeps10.dat");
    FILE *energy_file1 = fopen(energy_filename1, "w");

    printf("Starting Equilibration Phase...........\n");
    #pragma omp parallel  // Creates a parallel region once
    {
        // #pragma omp single
        // {
        //     int num_threads = omp_get_num_threads();
        //     printf("Number of threads activated for omp parallel: %d\n", num_threads);
        // }


    for (int S=0; S<EQUI_SWEEP; S++) {

        // #pragma omp single
        // {printf("[Sweep %d]running on thread %d............\n", 
        //                S, omp_get_thread_num());}


        // #pragma omp master
        // {
            #pragma omp for schedule(dynamic)
            for (int i = 0; i < NUM_REPLICAS; i++) {

            // printf("[Sweep %d] Replica %d running on thread %d: Energy = %f, Mean = %f\n", 
            //            S, i, omp_get_thread_num(),
            //            replicas[i]->H_e,
            //            replicas[i]->energy_stats.sample_mean);

            mc_steps(replicas[i]);

            update_running_stats(&replicas[i]->energy_stats, replicas[i]->H_e);

                }
        // }
        
        #pragma omp barrier

        // #pragma omp single
        // {printf("After %d sweeps:\n", S);}

        if (signal_was_caught()) {
                    printf("Signal caught. Exiting Replicas Equilibration Calculation...\n");
                    break;
                }


    }
    #pragma omp single
    {printf("Equilibration Phase Done...........\n");}
    

    // exit(0);

    // printf("All replicas equilibrated after %d sweeps.\n", S);
    
    #pragma omp single 
    { 
    printf("Starting Production Phase...........\n");
    printf("Initialise the Running Stats\n");

    for (int i = 0; i < NUM_REPLICAS; i++) {
            // mc_steps(replicas[i], p);
            // update_running_stats(&replicas[i]->energy_stats, replicas[i]->H_e);

            printf("Replica %d: Energy: %f, Mean: %f, Std Dev: %f, Std Err: %f\n", i, replicas[i]->H_e, replicas[i]->energy_stats.sample_mean, replicas[i]->energy_stats.sample_std_dev, replicas[i]->energy_stats.standard_error);
               
        }
    }

   

    // exit(0);

    for (int S= 0; S < TOT_SWEEP; S++) {
        // #pragma omp single 
        // {printf("Sweep %d, Running on Thread %d ...............\n", S, omp_get_thread_num());}

        #pragma omp for schedule(dynamic)
        for (int i = 0; i < NUM_REPLICAS; i++) 
        {
                
                mc_steps(replicas[i]); 
                fprintf(energy_file1, "%f\n", replicas[i]->H_e);

                // update_running_stats(&replicas[i]->energy_stats, replicas[i]->H_e);

            }
        #pragma omp barrier
        #pragma omp master 
        {
            if (S%100==0){

                    printf("After %d sweeps:\n", S);
                    
                    tries++;

                    // performReplicaSwaps(replicas, NUM_REPLICAS);

                    for (int i = 0; i < NUM_REPLICAS; i++) {

                        // double accept_rate = (double)replicas[i]->accepted_swaps/(double)replicas[i]->swap_tries;
                        
                        // printf("Replica %d: Energy: %f, Mean: %f, Std Dev: %f, Std Err: %f, Acc_Swap:%ld, Tried_Swap:%ld, Swap_rate: %f\n", i, replicas[i]->H_e, 
                        //                                                                             replicas[i]->energy_stats.sample_mean,
                        //                                                                             replicas[i]->energy_stats.sample_std_dev, 
                        //                                                                             replicas[i]->energy_stats.standard_error,
                        //                                                                             replicas[i]->accepted_swaps,
                        //                                                                             replicas[i]->swap_tries,
                        //                                                                             accept_rate);
                        // printf("\n");

                        double energy_diff = replicas[i]->energy_B - replicas[i]->energy_A;
                        energies[i] += energy_diff;

                        
                        }

                    // time_t now = time(NULL);
                    // char *timestamp = ctime(&now);
                    // timestamp[strlen(timestamp) - 1] = '\0'; // Remove the newline character
                    for (int i = 0; i < NUM_REPLICAS; i++) {
                       

                    //     double energy_diff = replicas[i]->energy_B - replicas[i]->energy_A;
                        // fprintf(energy_file, "%f\t%f\t%f\n", replicas[i]->mu, energy_diff, replicas[i]->H_e);
                        fprintf(energy_file, "%f\n", replicas[i]->H_e);


                    //     // printf("[%s] Replica Mu: %f: Energy: %f, Mean: %f, Std Dev: %f, Std Err: %f\n", timestamp, replicas[i]->mu, replicas[i]->H_e, 
                    //                                                                                 // replicas[i]->energy_stats.sample_mean,
                    //                                                                                 // replicas[i]->energy_stats.sample_std_dev, 
                    //                                                                                 // replicas[i]->energy_stats.standard_error);
                                                                                                 
                    
                    }
                    // // fprintf(energy_file, "\n");
                    // // printf("\n");

                    // fflush(energy_file);


                

            }
        }


        #pragma omp barrier

        if (signal_was_caught()) {
                    printf("Signal caught. Exiting Production Phase Calculation...\n");
                    break;

                        }
        

    }
    
    }
    
    
    
    // Free memory
    for (int i = 0; i < NUM_REPLICAS; i++) {
        freeReplica(replicas[i]);
    }


    gsl_rng_free (rg);

    if (tries == 0) {
        printf("No swaps were performed, cannot compute average energies.\n");
    } else {
        printf("Averaged Energies:\n");
        for (int i = 0; i < NUM_REPLICAS; i++) {
            double avg_energy = energies[i] / tries;
            printf("Replica %f: %f\n", Mus[i], avg_energy);
        }
    }

    gettimeofday(&tv_end, NULL);
    double elapsed_seconds = (tv_end.tv_sec - tv_start.tv_sec) 
                             + (tv_end.tv_usec - tv_start.tv_usec) / 1e6;
    printf("Elapsed wall-clock time: %f seconds\n", elapsed_seconds);


    fclose(energy_file);
    fclose(energy_file1);



}