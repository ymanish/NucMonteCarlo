// Header file

#include "dnaMC2.h"

// Includes
#include <signal.h>
#include <unistd.h>

// MC Sweep size
const int ENSSTEP = 147;

gsl_rng * r = NULL;

double sum_energy = 0.0;
double sum_energy_squared = 0.0;
double mean_energy, std_dev, std_err;
long int num_samples = 0;

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

static int signal_was_caught(void)
{
    if (sigint_received) fprintf(stderr, "SIGINT received!\n");
    if (sigterm_received) fprintf(stderr, "SIGTERM received!\n");
    if (sigquit_received) fprintf(stderr, "SIGQUIT received!\n");
    return (sigint_received || sigterm_received || sigquit_received);
}
// END OF SIGNAL HANDLING CODE



// Globally accessible sequence arrays
Base *Sequence;
Base *PullSequence;

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

void compute_bound_sites(int left, int right) {
    // Original boundbp array
    int boundbp_original[29] = {
        3, 7, 15, 18, 25, 30, 35, 39, 46, 50, 56, 60, 66, 70,
        77, 81, 87, 91, 97, 101, 108, 112, 117, 122, 129, 132, 140, 144, 5000
    };

    // Starting and ending indices for slicing
    int start = left * 2;                 // Multiply by 2 to account for pairs
    int end = (right * 2) + 1;            // Include all up to the right * 2 + 2
    if (start < 0 || start >= 29 || end <= 0 || end > 29 || start >= end) {
        fprintf(stderr, "Error: Invalid unwrapping parameters (left=%d, right=%d).\n", left, right);
        exit(1);
    }

    // Update the boundbp array
    for (int i = 0; i < 29; i++) {
        if (i >= start && i <= end) {
            boundbp[i] = boundbp_original[i]; // Keep original value in range
        } else {
            boundbp[i] = 5000; // Assign a large number for out-of-range values
        }
    }
}




// Return 1 if base is bound to base+1, -1 if to base-1, 0 if free.
// The function determines this based on the boundbp array defined above.
int isBound(int base) {    
    int base2 = base + 1;
    int n;
    
    // Copy base into the last (unused) element of boundbp
    boundbp[28] = base;
    
    // Loop over boundbp until we find the index at which
    // base is located
    for ( n = 0; boundbp[n] != base; n++ );
    
    // If base was not found until at the end of the array,
    // it's not present in the array.
    if ( n == 28 ) {
        // But boundbp only contains one bp for every set
        // of bound base pairs, so we also have to check
        // for base+1. We use the same procedure.
        boundbp[28] = base2;
        for ( n = 0; boundbp[n] != base2; n++ );
        
        // If the potential complement was also not found,
        // the base pair is not bound. Return 0.
        if ( n == 28 ) { return 0; }
        // If it was found, the base pair is bound to 
        // base+1. Return 1.
        else { return 1; }
    }
    
    // If the outer if-statement above triggered, the 
    // following line will never be reached. If we do
    // get here, it means base was found in boundbp.
    // Return -1.
    return -1;
    
}

// Convert Base number to Bond number (BaseToBond)
//
// This function takes the number of the base along a nucleosome and checks
// if that base is one of the ones listed in the boundbp array.
int BaseToBond(int base) {
    int n;
    for ( n = 0; n < 28; n++ ) {
        if ( boundbp[n] == base ) {
            return n;
        }
    }
    return -1;
}

// Procedure to make sure the given matrix is orthogonal.
// It does so by converting to axis-angle representation
// and back to a rotation matrix. If the rotation matrix 
// is not perfectly orthogonal, the axis and angle will
// also come out imperfectly. The procedure forcibly
// renormalizes the axis and assumes the angle is a good
// enough approximation. Should only be used on matrices
// that are very approximately orthogonal.


// The Role of a 3x3 Matrix in Simulations:

// A 3x3 matrix in such simulations often represents a rotation matrix. This matrix is used to rotate or transform vectors in three-dimensional space, which is essential for accurately modeling the spatial configuration of the nucleosome.
// Requirement for Orthogonality:

// For a rotation matrix, it's crucial that the matrix is orthogonal. An orthogonal matrix ensures that the length of vectors (i.e., their magnitude) and the angles between them are preserved after transformation. This is vital for maintaining the physical integrity of the nucleosome structure during rotation or transformation.


void OrthogonalizeR(double R[3][3]) {
    double a[3];
    int i;
    
    // Get rotation angle
    double rot_ang = acos(0.5*(R[0][0] + R[1][1] + R[2][2] - 1.0));
    
    // Get rotation axis. 
    a[0] = R[2][1] - R[1][2];
    a[1] = R[0][2] - R[2][0];
    a[2] = R[1][0] - R[0][1];
    
    // Normalize axis.
    double t = 1.0/sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
    for ( i = 0; i < 3; i++ ) {
        a[i] *= t;
    }   
    
    // Rebuild rotation matrix. Store in same array.
    RfromAA(rot_ang, a, R);
}


// Signal handling function
void signal_handler(int signal_num) {
    printf("Interrupt signal is (%d).\n", signal_num);
    exit(signal_num);  
}

void printDNAState(DNAState state){
        printf("Position: (%f, %f, %f)\n", state.pos[0], state.pos[1], state.pos[2]);

        printf("Rotation Matrix:\n");
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                printf("%f ", state.R[i][j]);
            }
            printf("\n");
        }
}



void saveStateToPDB(const char *filename, DNAState state[N]) {
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        fprintf(stderr, "Error: Could not open file %s for writing.\n", filename);
        return;
    }

    fprintf(file, "REMARK  Generated by Nucleosome MMC\n");

    // Write each base pair as a single "atom" in the PDB file
    for (int i = 0; i < N; i++) {
        fprintf(file,
                "ATOM  %5d  P   DNA A%4d    %8.3f%8.3f%8.3f  1.00  0.00\n",
                i + 1,                              // Atom serial number
                i + 1,                              // Residue number
                state[i].pos[0],                    // X coordinate
                state[i].pos[1],                    // Y coordinate
                state[i].pos[2]);                   // Z coordinate
    }

    fprintf(file, "END\n");
    fclose(file);
}


// Save DNAState to an .xyz file
void saveStateToXYZ(const char *filename, DNAState state[N]) {
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        fprintf(stderr, "Error: Could not open file %s for writing.\n", filename);
        return;
    }

    // Write the number of atoms
    fprintf(file, "%d\n", N);

    // Write a blank line for the molecule name (optional)
    fprintf(file, "Nucleosome State\n");

    // Write atom positions in the format: "AtomName X Y Z"
    for (int i = 0; i < N; i++) {
        fprintf(file, "P  %.3f  %.3f  %.3f\n",  // Using "P" for phosphate
                state[i].pos[0],               // X coordinate
                state[i].pos[1],               // Y coordinate
                state[i].pos[2]);              // Z coordinate
    }

    fclose(file);
}


void cranklength(int m , int *idA, int *idB) {

    int mdfrom = boundbp[m]+1;
    int mdto = boundbp[m+1]-1;
    int hinge_dist = 0;
    *idA = 0;
    *idB = 0;

    printf("Insid M %d\n", m);

    // printf("#############################");
    // printf("CASE mdfrom %d , mdto %d\n", mdfrom, mdto);

    if (mdto - mdfrom >= 3) {
        // printf("Big enough\n");

        if (mdto - mdfrom == 3) {
            *idA = mdfrom;
            *idB = mdto;
            // printf("Special hinge dist=3  idA %d , idB %d\n", *idA, *idB);

        }
        else {
            // printf("Large Big enough\n");
            // printf("diff1 %d\n", (mdto - mdfrom)-2);
            // printf("mdfrom %d\n", mdfrom);
            // printf("random %d\n", d);
            *idA = mdfrom + gsl_rng_uniform_int(r, (mdto - mdfrom)-2);
            // printf("idA %d\n", *idA);
            

            hinge_dist = 3 + gsl_rng_uniform_int(r, (mdto - *idA -2));
            // printf("diff2 %d\n", (mdto - *idA -2));
            // printf("hinge_dist %d\n", hinge_dist);

            *idB = *idA + hinge_dist;
            if (*idB > mdto) {
                // printf("idB > mdto\n");
                *idB = mdto;
            }

            // printf("Hinge distance: %d, idA %d , idB %d\n", hinge_dist, *idA, *idB);
            }
        
    }}




// Main function, starting point of the program. The arguments
// encode the number and values of the command line arguments
// passed when calling the program.
int main(int argc, char *argv[])
{
    // Start signal handling
    // setup_signal_handler(SIGINT, handle_sigint);
    // setup_signal_handler(SIGTERM, handle_sigterm);
    // setup_signal_handler(SIGQUIT, handle_sigquit);
    signal(SIGINT, signal_handler);
    
    // Loop counters
    int i, n, n1;
    
    int idA;
    int idB;

    // Auxiliary
    int accept, tries;
    
    // Temporary base pair identity storage
    Base stor;
    
    // Step size along long sequence
    int dec = 1;
    
    // MC loop counter and length
    long int c;
    long int stop = 1e4;
    
    // Inverse temperature
    double beta = 1;
    
    // Equipartition energy, currently not used
    // double Eth = 3.0*(147.0-28.0)/beta;
    
    // DNAState arrays
    DNAState state[N];
    DNAState start[N];
    DNAState trial[N];
    
    // Midplane arrays
    DNAState midstate[N-1];
    DNAState midplanes[28];
    
    // BPS energy bookkeeping arrays. These arrays record the
    // local elastic energies in all base pair steps in the system
    // for the current state of the simulation, as well as for the
    // trial state under consideration during a move attempt.
    // The reason for these arrays' existence is that in order to
    // calculate the energy difference of a trial move with the current
    // state we need the energy of both the trial and current states.
    // However, we can reduce the number of calculations required by
    // realizing the following:
    // 1) When we accept a move, the new "current" energy becomes the
    //    trial energy. We can save calculations in the next move by 
    //    remembering this energy.
    // 2) We perform many localized moves that only affect the energy
    //    of a small part of the system. By breaking up our recollection
    //    of the energy down to the BPS level, we can recalculate and
    //    update only the base pair steps affected by a move.
    // The above two optimizations are addressed by introducing these
    // two bookkeeping arrays. The updating of the Trial array is done
    // within the various Move functions, copying over of the appropriate
    // values is done below in the MC loop.
    double BPSEnergies[N-1], BPSEnergiesTrial[N-1];
    
    // Stepsize for random MC moves
    double stepsize_t = 1.0/(sqrt(beta)*15.0); 
    double stepsize_a = PI/(sqrt(beta)*150.0);
    // double crankstepsize_a = 1.0/(sqrt(beta));
    
    // Variables for the MC loop:
    // e  - energy of current state
    // d  - random floating point container
    // w  - Boltzmann weight container
    // dE - energy difference for trial move
    // s  - random floating point for move type selection
    double e, d, w, dE, s, e_af, e_bf;
    
    // Mutation move probability
    // double Pmut = 0.2;
    double Pmut = 0.0;

    
    // Variables related to bound base pairs
    int bb, bplo, bphi;
    
    // Random seed for the random number generator.
    // It is set here to an arbitrary number. If a
    // seed value is provided on the command line,
    // this value is overwritten. Not providing a
    // value will lead to this same value always
    // being used, with reproducible results.
    int seed = 10089;
    
    // Start position and number of base pairs
    // to indicate region along long sequence
    // to analyse
    int startpos, nstep;

    int base_sel_left, base_sel_right;
    
    // Auxiliary variables for reading in sequences
    int numCharacters = 0, cnt = 0;
    
    // Read in sequence start position and number of steps 
    // to take from the command line. See README for full
    // definition of the command line arguments. See comments
    // describing the function of the two main nested loops
    // below for a description of how these variables are used.
    if (argc > 1) {
        startpos = atoi(argv[1]);
        nstep = atoi(argv[2]);

        if (argc > 3) {
            base_sel_left = atoi(argv[3]);
            base_sel_right = atoi(argv[4]);
        }
    }
    else { // Default is to start at the beginning and take only the first position
        startpos = 0;
        nstep = 1;
        base_sel_left = 0;
        base_sel_right = 13;

    }
    
    compute_bound_sites(base_sel_left, base_sel_right);

    // Debug: Print the modified boundbp array
    printf("Modified boundbp array:\n");
    for (int i = 0; i < 29; i++) {
        printf("%d ", boundbp[i]);
    }
    printf("\n");

    // exit(0);
    // Set the nucleosome-bound sequence to be the 147-bp sequence starting at 
    // the starting position
    Sequence += startpos;
    
    // Initialization routines
    init(seed);

    // exit(0);

    
    fillK();
    
    // Read in the sequence once to determine the length
    fprintf(stderr, "Reading sequence length.\n");
    FILE* seqfile = fopen("./State/SequenceA.seq", "r");
    char nextChar = getc(seqfile);
    while (nextChar != EOF) {
        if ( isspace(nextChar) ) { nextChar = getc(seqfile); }
        else {
            numCharacters++;
            nextChar = getc(seqfile);
        }
    }
    fclose(seqfile);
    
    // Allocate memory to hold the sequence
    fprintf(stderr, "Allocating space for sequence.\n");
    PullSequence = (Base*) malloc(numCharacters*sizeof(int));
    Sequence = PullSequence;
    
    // Read in the sequence again and store in memory
    fprintf(stderr, "Reading sequence.\n");
    seqfile = fopen("./State/SequenceA.seq", "r");
    nextChar = getc(seqfile);
    while (nextChar != EOF) {
        if ( isspace(nextChar) ) { nextChar = getc(seqfile); }
        else {
            Sequence[cnt] = CharToBase(nextChar);
            cnt++;
            nextChar = getc(seqfile);
        }
    }
    fclose(seqfile);
    fprintf(stderr, "Sequence read.\n");
    // printf("%d", numCharacters);
    // printf("%zu\n", sizeof(Sequence));
    // exit(0);
    // Read in nucleosome configuration file
    fprintf(stderr, "Reading state.\n");
    FILE* statefile = fopen("State/Nucleosome.state", "r");



    ReadDNAState(statefile, state);


    // exit(0);

    fclose(statefile);
    fprintf(stderr, "State read.\n");
    
    
    printDNAState(state[0]);  // Print the contents of myState
   
    // Calculate midplanes for the state just read in
    midplanestate(state, midstate);
    
    // Store the 28 relevant midframes in a new array
    for ( i = 0; i < N; i++ ) {
        if ( BaseToBond(i) >= 0 ) {
            //fprintf(stderr, "%d\n", i);
            midplanes[BaseToBond(i)] = midstate[i-1];
        }
    }
    
    // Set the starting configuration for the MC loops
    // to the state just read in
    CopyDNAState(state, start, 0, N-1);
    
    // Starting energy with the given configuration
    // and sequence.
    e = Energy(state, Sequence);
    fprintf(stderr, "AVG Energy: %e\n", e);
    
    // Done initializing, let's start
    fprintf(stderr, "Starting MC.\n");
    

    // exit(0);

    FILE *energy_file = fopen("energies.dat", "w");
    if (energy_file == NULL) {
        fprintf(stderr, "Error: Could not open energies.dat for writing.\n");
        exit(1);
    }


    // The program has two nested loops:
    //
    // 1) Loop with counter i. Starting at position startpos along a
    // sequence of length at least 147 + startpos + nstep, it
    // loops over all positions (startpos + n*dec) < startpos + nstep
    // along the given sequence and runs MC for the 147-bp subsequence
    // at each position
    //.
    // 2) The actual MC loop with counter c. It runs a number of sweeps
    // defined by the variable stop. In each sweep, it makes ENSSTEP moves.
    // Only accepted moves are counted.
    
    for ( i = startpos; i < startpos+nstep; i+=dec ) {
        // Initialize variables
        

        // Always start in the same state
        CopyDNAState(start, state, 0, N-1);
        


         // Start energy tracking
        for ( n = 0; n < N; n++ ) {
            // For safety, orthogonalize the orientations before starting
            OrthogonalizeR(state[n].R);
        }

        // Clean trial state
        CopyDNAState(state, trial, 0, N-1);
              

        //// Here the state are not orthogonalized 
        e = Energy(state, Sequence);
        e_af = e;
        
        for ( n = 0; n < N-1; n++ ) {
            // For safety, orthogonalize the orientations before starting
            // OrthogonalizeR(state[n].R);
            
            // Set up bookkeeping arrays
            BPSEnergies[n] = BPSEnergy(state[n], state[n+1], Sequence[n], Sequence[n+1]);
            BPSEnergiesTrial[n] = BPSEnergies[n];
        }
        
        // Reset number of attempted MC moves
        tries = 0;
        
        // Start MC loop
        // We make a number of sweeps defined by stop
        // Each sweep consists of ENSSTEP moves
        for ( c = 0; c < (ENSSTEP*stop); c++ ) {
            
            // Loop until accepted move is found
            accept = 0;
            while ( accept == 0 ) {
                // Increase attempt number
                tries += 1;
                
                // Select random base pair
                n1 = gsl_rng_uniform_int(r, N);     
                
                // Check if base pair is bound
                bb = isBound(n1);                   
                
                // If the base pair is bound, bb is either +1 or -1,
                // indicating whether it is bound to its left or right
                // neighbour. The following ternary expressions lead to
            // the following results, depending on this bound state:
                //
                // Case 1: n1 is not bound -> bplo = bphi = n1
                // Case 2: n1 is bound to the left neighbour -> bplo = n1-1, bphi = n1
                // Case 3: n1 is bound to the right neighbour -> bplo = n1, bphi = n1+1
                //
                // The bphi and bplo variables thus encode which base pairs
                // need to be moved, and thus which energies neeed to be
                // recalculated.
                //
                // Explanation of the ternary expressions:
                // 
                // bplo : if n1 < n1 + bb, then the base pair is bound to its right neighbour,
                // and bplo needs to be n1. If not, then there are two possibilities:
                // 1) bb = 0  : base pair not bound, bplo should be equal to n1, but we can safely add bb as it is 0.
                // 2) bb = -1 : base pair bound on the left, bplo needs to be n1 - 1 = n1 + bb.
                // Thus, in both cases bplo = n1 + bb.
                bplo = ( n1 < n1+bb ? n1 : n1+bb ); 
                // bphi : analogous to the above. If n1 < n1 + bb, we are bound on the right, so we
                // need to add bb to n1. In the other two cases, bphi = n1.
                bphi = ( n1 < n1+bb ? n1+bb : n1 );
                
                // Select random move type
                s = gsl_rng_uniform(r);     
                
                // Store the sequence identity of the selected base pair
                // (only needed in case of mutation).    
                stor = Sequence[n1];
                
                // Select mutation with probability Pmut
                if ( s < Pmut ) {
                    dE = Mutation(state, Sequence, n1, &stor, BPSEnergies, BPSEnergiesTrial);
                }
                // Else, we do a configurational move. If the base pair 
                // is bound, keep midplane fixed.
                else if ( bb ) { 
                    
                    // printf("Start of SBPMOVE FMP move\n");
                    // printf(e_bf);
                    // e_bf = Energy(state, Sequence);


                    dE = SBPMoveFMP(bplo, bphi, stepsize_a, stepsize_t, state, trial, BPSEnergies, BPSEnergiesTrial, midplanes);
                    
                }
                // Otherwise, move freely.
                else {      
                    dE = SBPMove(n1, 2.0*stepsize_a, 2.0*stepsize_t, state, trial, BPSEnergies, BPSEnergiesTrial); 
                }
                
                // Acceptance conditions.
                // If the move lowers the energy, we always accept it.
                if ( dE < 0 ) { 
                    // Setting accept = 1 will end the while loop after this iteration
                    accept = 1;
                    // printf("Accepted move\n");
                    
                    // Copy the trial state to the current state
                    CopyDNAState(trial, state, bplo, bphi);
                    
                    // Copy the trial energies as well
                    CopyBPSEnergies(BPSEnergiesTrial, BPSEnergies, bplo-1, bphi);
                    
                    // Copy the mutation (does nothing if we didn't mutate)
                    Sequence[n1] = stor;
                    
                    // Update the energy
                    e += dE;
                }
                // If the move doesn't lower the energy, we accept it with
                // probability exp(-beta*dE).
                else { 
                    // Random float between 0 and 1.
                    d = gsl_rng_uniform(r);
                    
                    // Boltzmann weight of the move.
                    w = exp(-beta*dE); 
                    
                    // With probability given by w, accept the move anyway.
                    if ( d < w ) {
                        // Same acceptance operations as above.
                        accept = 1;

                        // printf("Accepted move with boltzmann\n");

                        CopyDNAState(trial, state, bplo, bphi);
                        CopyBPSEnergies(BPSEnergiesTrial, BPSEnergies, bplo-1, bphi);
                        e += dE;
                        Sequence[n1] = stor;
                    }
                    // With probability 1-w, reject the move.
                    else {

                        // printf("Rejected move......\n");
                        // Clean up trial state
                        CopyDNAState(state, trial, bplo, bphi);
                        
                        // Clean up bookkeeping
                        CopyBPSEnergies(BPSEnergies, BPSEnergiesTrial, bplo-1, bphi);
                    }
                }

                double norm_diff; 
                double diff[3] = {0.0, 0.0, 0.0};
                    
                for (int m = 0; m<27; m++) {
                    
                    printf("Base Pair %d\n", m);
                    printf("Boundbp %d\n", boundbp[m]);
                    printf("Boundbp %d\n", boundbp[m+1]);

                    cranklength(m , &idA, &idB);

                    if (idA == 0 || idB == 0) {continue;}
                                                    
                    else {
                        // printf("idA %d , idB %d\n", idA, idB);
                        /// Run the crankshaft move

                        // printf("Makig MOVEEEEEEEEEEEEE\n");

                        // for (int i = 0; i <N; i++) {
                        //     check_orthogonal(state[i].R);
                        // }
                        

                        e_bf = Energy(state, Sequence);
                        
                        // printf("Energy before move %f\n", e_bf);


                        // printf("Energy per Base Pair before move.......\n");
                        // // printf("dE %f\n", dE);
                        // for (int i = idA-1; i <= idB; i++) {
                        //     printf("BPSEnergiesTrial[%d] = %f\n", i, BPSEnergiesTrial[i]);
                        //     printf("BPSEnergies[%d] = %f\n", i, BPSEnergies[i]);
                        // }


                        
                     
                        // for (int i = idA-2; i <= idB+2; i++) {
                        //     printf("Base Number %d\n", i);
                        //     PrintBPState_terminal(state[i]);

                        // }


                        // for (int k = idA; k <= idB; k++) {for (i = 0; i < 3; i++) {diff[i] = state[k+1].pos[i] - state[k].pos[i];}
                        //     printf("ID %d  norm_diff %f\n", k, vec_norm(diff));
                        //     printf("Relative Angle %f\n", realtive_angle(state[k], state[k+1]));
                        // }

                    


                        dE = CrankShaftMove(idA, idB, stepsize_a, beta, state, trial, BPSEnergies, BPSEnergiesTrial);

                        // printf("Move has been made...................\n");

                        // for (int i = idA-1; i <= idB-1; i++) {
                        //     printf("BPSEnergiesTrial[%d] = %f\n", i, BPSEnergiesTrial[i]);
                        // }





                        CopyDNAState(trial, state, idA, idB);
                        // CopyBPSEnergies(BPSEnergiesTrial, BPSEnergies, idA-1, idB-1);
                        BPSEnergies[idA-1] = BPSEnergiesTrial[idA-1];
                        BPSEnergies[idB-1] = BPSEnergiesTrial[idB-1];

                        // for (int i = idA-2; i <= idB+2; i++) {
                        //     printf("Base Number %d\n", i);
                        //     PrintBPState_terminal(state[i]);
                        // }

                        // for (int i = idA; i <= idB; i++) {
                        //     PrintBPState_terminal(state[i]);
                        // }

                        // for (int k = idA; k <= idB; k++) {for (i = 0; i < 3; i++) {diff[i] = state[k+1].pos[i] - state[k].pos[i];}
                        //     printf("ID %d  norm_diff %f\n", k, vec_norm(diff));
                        //     printf("Relative Angle %f\n", realtive_angle(state[k], state[k+1]));
                        // }




                        // printf("Energy per Base Pair after move.......\n");
                        // // printf("dE %f\n", dE);
                        // for (int i = idA-1; i <= idB; i++) {
                        //     printf("BPSEnergiesTrial[%d] = %f\n", i, BPSEnergiesTrial[i]);
                        //     printf("BPSEnergies[%d] = %f\n", i, BPSEnergies[i]);
                        // }



                        // printf("BPSEnergies after move\n");
                        //     for (int i = idA-1; i <= idB; i++) {
                        //     // OrthogonalizeR(state[i].R);
                        //     BPSEnergies[i] = BPSEnergy(state[i], state[i+1], Sequence[i], Sequence[i+1]);
                        //     printf("BPSEnergies[%d] = %f\n", i, BPSEnergy(state[i], state[i+1], Sequence[i], Sequence[i+1]));
                        //     printf("BPSEnergies[%d] = %f\n", i, BPSEnergies[i]);
                        // }
                        // for ( n = 0; n < N; n++ ) {
                        // // For safety, orthogonalize the orientations 
                        // OrthogonalizeR(state[n].R);
                        // }
                    
                        
                        // for (int i = 0; i <N; i++) {
                        //     check_orthogonal(state[i].R);
                        // }
                        
                        e_af = Energy(state, Sequence);

                        // printf("Energy after move..... %f\n", e_af);
                    
                    
                        printf("Step number %ld\n", c);
                        printf("Energy difference..... %.15f\n", e_af-e_bf);
                        printf("Energy difference dE.. %.15f\n", dE);



                        // printf("\n");
                        // printf("\n");
                        // printf("\n");
                        // printf("\n");
                     
                        e += dE;
                        // if (fabs(e_af - e_bf - dE) > 1e-8) {
                        //     printf("Energy difference not equal to dE\n");
                        //     exit(0);
                        // }
                        // exit(0);
                    }
                }
                // exit(0);

                if (tries%100==0){
                    fprintf(energy_file, "%ld\t%.6f\n", c, e);
                    
                }
                
            }
            
            // We usually want to output something during the simulation.
            // Here are some examples. The current statement will print
            // the sweep number, the tracked energy, and a freshly calculated
            // energy. The latter two should be identical, otherwise there
            // is a leak somewhere.
            // The lines that are commented out would instead print the
            // full DNA configuration or the sequence.
            if ( c % (1*ENSSTEP*1000) == 0 ) {
                fprintf(stdout, "%ld\t%e\t%e\n", c/ENSSTEP, e, Energy(state, Sequence));
                
                // for ( n = 0; n < N; n++ ) {
                //     // For safety, orthogonalize the orientations 
                //     OrthogonalizeR(state[n].R);
                // }
                

                // PrintDNAState(stdout, state);
                // PrintSequence(stdout, Sequence);

                char filename[100];
                // snprintf(filename, sizeof(filename), "./Trajectories/Nucleosome_%ld.pdb", c / (ENSSTEP * 1000));
                // saveStateToPDB(filename, state);  // Save the state to a PDB file

                snprintf(filename, sizeof(filename), "./Trajectories/Frame_%ld.xyz", c / (ENSSTEP * 1000));
                saveStateToXYZ(filename, state);  // Save the state to a PDB file
                fprintf(stderr, "Saved state to %s\n", filename);
                fflush(stdout);
            }
            
            // If any signal was caught, break out of the loop
            if (signal_was_caught(  )) break;



            // sum_energy += e;
            // sum_energy_squared += e * e;
            // num_samples++;

            // if (c > 100){
            //     // Calculate mean, standard deviation, and standard error
            //     mean_energy = sum_energy / num_samples;
            //     std_dev = sqrt((sum_energy_squared / num_samples) - (mean_energy * mean_energy));
            //     std_err = std_dev / sqrt(num_samples);
            //     fprintf(stderr, "Mean: %.6f, Standard Deviation: %.6f, Standard Error: %.6f\n", mean_energy, std_dev, std_err);
            //     // Check if standard error is below the threshold
            //     if (std_err < 0.01) {
            //         fprintf(stderr, "Stopping simulation: Standard error below threshold (%.6f).\n", std_err);
            //         break;  // Exit the loop
            //     }
            // }
           
        }
        
        combineXYZFiles();
        // Close the file
        fclose(energy_file);

        // If any signal was caught, break out of the loop
        if (signal_was_caught(  )) break;
    }
    
    // Clean up the random number generator.
    gsl_rng_free (r);
    
    // Cleanly exit, returning exit status 0 (success).
    return 0;
}
