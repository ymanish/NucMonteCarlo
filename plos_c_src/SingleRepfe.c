// Header file

#include "dnaMC2.h"
#include <sys/time.h>
#include <signal.h>
#include <unistd.h>
#include <stdbool.h>

// MC Sweep size
const int ENSSTEP = 147;

gsl_rng * rg = NULL;
// gsl_rng * r2 = NULL;

// double betas[NUM_REPLICAS] = {1.0,0.9,0.8,0.7};

// double Mus[NUM_REPLICAS] = {0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0};
// double Mus[NUM_REPLICAS] = {0.5};

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
    
    int i, n, n1, l, id_S, id_E, ID;
    int seed = 10089; 
    int startpos, nstep;
    int base_sel_left, base_sel_right;
    double beta, gmu;
    double e_A, e_B, H_e;
    long int TOT_SWEEP = 1e5;
    long int EQUI_SWEEP = 1e4;
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

    const char *seqA_str;
    const char *seqB_str;

    int SINGLE_REP_INDEX = 0;

    if (argc > 1) {
        startpos = atoi(argv[1]);
        nstep = atoi(argv[2]);

        if (argc > 3) {
            base_sel_left = atoi(argv[3]);
            base_sel_right = atoi(argv[4]);
        }
        if (argc > 4) {
            gmu = atof(argv[5]);
            seqA_str = argv[6];
            seqB_str = argv[7];
            ID = atoi(argv[8]);

        }
    }
    else { // Default is to start at the beginning and take only the first position
        startpos = 0;
        nstep = 1;
        base_sel_left = 0;
        base_sel_right = 13;
        gmu = 1;
        seqA_str = "CCGCTCAATTGGTCGTAGACAGCTCGCACGTACGCCTGGAGAATCCCGGTGCCGAGGGCTGTCCCCCGCGTTTCCCTAGTCTTAACTCCAGGCACGTGTTAGCACCGCTTAAACCAGATATATACATCCTGTCGCCAAGGGGATTAC"; // Random Seq
        seqB_str = "CTGGAGAATCCCGGTGCCGAGGCCGCTCAATTGGTCGTAGACAGCTCTAGCACCGCTTAAACGCACGTACGCGCTGTCCCCCGCGTTTTAACCGCCAAGGGGATTACTCCCTAGTCTCCAGGCACGTGTCAGATATATACATCCTGT"; //601 Seq
        ID = 1;
    }
    
    // Get sequences from command line
   

    // Determine lengths
    size_t lenA = strlen(seqA_str);
    size_t lenB = strlen(seqB_str);

    
    Sequence += startpos;
    Sequence2  += startpos;

    // Initialization routines
    init(seed);
    // init(seed+1);

    
    fillK();


    // Read in sequences from SequenceA and SequenceB files
    // Allocate memory for sequences
    Sequence = (Base*) malloc(lenA * sizeof(Base));
    if (Sequence == NULL) {
        fprintf(stderr, "Memory allocation error for SequenceA.\n");
        exit(1);
    }
    Sequence2 = (Base*) malloc(lenB * sizeof(Base));
    if (Sequence2 == NULL) {
        fprintf(stderr, "Memory allocation error for SequenceB.\n");
        exit(1);
    }

    // Convert and store the sequences
    for (size_t i = 0; i < lenA; i++) {
        if (!isspace(seqA_str[i])) {
            Sequence[i] = CharToBase(seqA_str[i]);
        }
    }
    for (size_t i = 0; i < lenB; i++) {
        if (!isspace(seqB_str[i])) {
            Sequence2[i] = CharToBase(seqB_str[i]);
        }
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
        replicas[i]->mu = gmu;

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

     
         // Allocate a unique RNG for each replica.
        replicas[i]->rn = gsl_rng_alloc(gsl_rng_mt19937);
        gsl_rng_set(replicas[i]->rn, seed + i);

        compute_bound_sites(base_sel_left, base_sel_right, replicas[i]->boundbp);


    }


    // char energy_filename[200];
    // snprintf(energy_filename, sizeof(energy_filename), "./TI_output/Auto_Corr_SingleRep_S0_ssize10.dat");
    // FILE *energy_file = fopen(energy_filename, "w");
    // if (energy_file == NULL) {
    //     fprintf(stderr, "Error: Could not open energies.dat for writing.\n");
    //     exit(1);
    // }

    char energy_filename[200];
    snprintf(energy_filename, sizeof(energy_filename), "./TI_output/OneRep_OneMu/%d.dat", ID);
    FILE *energy_file = fopen(energy_filename, "w");


    // char energy_filename1[200];
    // snprintf(energy_filename1, sizeof(energy_filename1), "./TI_output/Auto_Corr_SingleRep_S100_ssize10.dat");
    // FILE *energy_file1 = fopen(energy_filename1, "w");



    printf("Starting Equilibration Phase...........\n");


    for (int S=0; S<EQUI_SWEEP; S++) {

        // mc_steps_shuffle(replicas[SINGLE_REP_INDEX]);
        mc_steps(replicas[SINGLE_REP_INDEX]); 

        // update_running_stats(&replicas[SINGLE_REP_INDEX]->energy_stats, replicas[SINGLE_REP_INDEX]->H_e);

    


        if (signal_was_caught()) {
                    printf("Signal caught. Exiting Replicas Equilibration Calculation...\n");
                    break;
                }


    }
    printf("Equilibration Phase Done...........\n");
    
    printf("Starting Production Phase...........\n");
    printf("Initialise the Running Stats\n");

    printf("Replica %d: Energy: %f, Mean: %f, Std Dev: %f, Std Err: %f\n", i, replicas[SINGLE_REP_INDEX]->H_e, replicas[SINGLE_REP_INDEX]->energy_stats.sample_mean, replicas[SINGLE_REP_INDEX]->energy_stats.sample_std_dev, replicas[SINGLE_REP_INDEX]->energy_stats.standard_error);
               
     

    for (int S= 0; S < TOT_SWEEP; S++) {
       
      
        // mc_steps_shuffle(replicas[SINGLE_REP_INDEX]); 
        mc_steps(replicas[SINGLE_REP_INDEX]); 

        fprintf(energy_file, "%f\n", replicas[SINGLE_REP_INDEX]->H_e);

        if (S%100==0){

            printf("After %d sweeps:\n", S);
            
            tries++;
            double energy_diff = replicas[SINGLE_REP_INDEX]->energy_B - replicas[SINGLE_REP_INDEX]->energy_A;
            energies[SINGLE_REP_INDEX] += energy_diff;
            
            update_running_stats(&replicas[SINGLE_REP_INDEX]->energy_stats, replicas[SINGLE_REP_INDEX]->H_e);

            // fprintf(energy_file1, "%f\n", replicas[SINGLE_REP_INDEX]->H_e);                                                          
            
        }
            

        
    


        if (signal_was_caught()) {
                    printf("Signal caught. Exiting Production Phase Calculation...\n");
                    break;

                        }
    
    }
    
    
    
    // Free memory
    for (int i = 0; i < NUM_REPLICAS; i++) {
        freeReplica(replicas[i]);
    }


    gsl_rng_free (rg);

    printf("Averaged Energies:\n");
    for (int i = 0; i < NUM_REPLICAS; i++) {
        double avg_energy = energies[i] / tries;
        printf("Mu Value: %f, Average Energy %f\n", replicas[i]->mu, avg_energy);

        // fprintf(energy_file, "%d,%f,%f\n", ID, replicas[i]->mu, avg_energy);


    }


    printf("ID: %d, Replica Mu: %f, Energy: %f, Mean: %f, Std Dev: %f, Std Err: %f\n", ID, replicas[SINGLE_REP_INDEX]->mu,
                                                                                replicas[SINGLE_REP_INDEX]->H_e, 
                                                                                replicas[SINGLE_REP_INDEX]->energy_stats.sample_mean,
                                                                                replicas[SINGLE_REP_INDEX]->energy_stats.sample_std_dev, 
                                                                                replicas[SINGLE_REP_INDEX]->energy_stats.standard_error);
                            


    gettimeofday(&tv_end, NULL);
    double elapsed_seconds = (tv_end.tv_sec - tv_start.tv_sec) 
                             + (tv_end.tv_usec - tv_start.tv_usec) / 1e6;
    printf("Elapsed wall-clock time: %f seconds\n", elapsed_seconds);


    fclose(energy_file);
    // fclose(energy_file1);



}