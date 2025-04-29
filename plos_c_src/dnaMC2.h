// Include all the libraries we need across different files
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <string.h>
#include <ctype.h>
#include <gsl/gsl_blas.h>
// #include <cblas.h>


// Defining some constants
#define PI 3.14159265358979

// System size needs to be hardcoded here
#define N 147

// The current code does not do tension, but
// I'll leave this here for the future.
#define f1pN 0.02414  // Pulling force of 1pN/(kT_r) * (m/Angstrom);
                      // 1pN of force, defined to give E/(kT) for room temperature if distance units are in Angstrom
                      

#define NUM_REPLICAS 1

// Define move types
#define MOVE_LOCAL       0
#define MOVE_LED         1
#define MOVE_RED         2
#define MOVE_CRANKSHAFT  3
#define MOVE_CLUSTER     4

// Global random number generator
extern gsl_rng * rg;  
// extern gsl_rng * r2;

// DNAState structure
typedef struct {
    double pos[3];
    double R[3][3];
} DNAState;

// Base variable type
typedef enum Bases {
	A, T, C, G, X
} Base;

typedef struct {
    long int sample;
    double value_sum;
    double sample_mean;
    double sum_sq_diff;
    double sample_variance;
    double sample_std_dev;
    double standard_error;

} RunningStats;

typedef struct {
    double beta;
    DNAState *state;
    DNAState *trial;
    DNAState *midstate;
    DNAState *midplanes;
    double energy_A;
    double energy_B;
    double H_e;
    double *BPSEnergies;
    double *BPSEnergies2;
    double *BPSEnergiesTrial;
    double *BPSEnergiesTrial2;
    long int accepted_swaps; 
    long int swap_tries;
    int equilibrated;
    double mu;
    gsl_rng *rn;
    RunningStats energy_stats;
    int *boundbp;
} Replica;




// Some variables that are defined in the C files but
// need to be available across files.
extern double Eq[4][4][6];                     
extern double Kf[4][4][6][6];
extern Base *Sequence;
extern char *BaseString[5];
extern int boundbp[29];

// For Thermodynamic Integration
extern Base *Sequence2;


// Function prototypes
int signal_was_caught(void);
void fillK();

void performReplicaSwaps(Replica **replicas, int num_replicas);

// Move functions
double SBPMove(int base, double stepsize_t, double stepsize_a, DNAState currentState[N], DNAState newState[N], double BPSEnergies[N-1], double BPSEnergiesTrial[N-1]);
double SBPMoveFMP(int base, int complement, double stepsize_t, double stepsize_a, DNAState currentState[N], DNAState newState[N], double BPSEnergies[N-1], double BPSEnergiesTrial[N-1], DNAState midplanes[28]);
double Mutation(DNAState state[N], Base* seq, int a, Base* b, double BPSEnergies[N-1], double BPSEnergiesTrial[N-1]);
double CrankShaftMove(int idA, int idB, double stepsize_a, double b, DNAState currentState[N], DNAState newState[N], double BPSEnergies[N-1], double BPSEnergiesTrial[N-1]);

// Move functions
double SBPMove2(int base, double stepsize_t, double stepsize_a, DNAState currentState[N], DNAState newState[N], double BPSEnergies[N-1], double BPSEnergiesTrial[N-1], double BPSEnergies2[N-1], double BPSEnergiesTrial2[N-1], double p, double *de_a, double *de_b, gsl_rng *r);
double SBPMoveFMP2(int base, int complement, double stepsize_t, double stepsize_a, DNAState currentState[N], DNAState newState[N], 
                    double BPSEnergies[N-1], double BPSEnergiesTrial[N-1], double BPSEnergies2[N-1], double BPSEnergiesTrial2[N-1], 
                    DNAState midplanes[28], double p, double *de_a, double *de_b, gsl_rng *r, int boundbp[29]);


double CrankShaftMove2(int id_S, int id_E, double stepsize_a, double beta, DNAState currentState[N], DNAState newState[N], double BPSEnergies[N-1], 
                    double BPSEnergiesTrial[N-1], double BPSEnergies2[N-1], double BPSEnergiesTrial2[N-1], double p, double *de_a, double *de_b, gsl_rng *r);


double ClusterTrans2(int id_S, int id_E, double stepsize_t, double beta, DNAState currentState[N], DNAState newState[N], 
                        double BPSEnergies[N-1], double BPSEnergiesTrial[N-1],
                         double BPSEnergies2[N-1], double BPSEnergiesTrial2[N-1], 
                         double p, double *de_a, double *de_b, gsl_rng *r);

void perform_sbp_move(int n1,
                        double stepsize_a,
                        double stepsize_t,
                        double beta,
                        DNAState state[N],
                        DNAState trial[N],
                        double BPSEnergies[N-1],
                        double BPSEnergiesTrial[N-1],
                        double BPSEnergies2[N-1],
                        double BPSEnergiesTrial2[N-1],
                        DNAState midplanes[28],
                        double p,
                        double *e_A,
                        double *e_B,
                        double *H_e,
                        long int *sbp_tries,
                        long int *midbp_tries,
                        long int *sbp_accept,
                        long int *midbp_accept, 
                        gsl_rng *r, int boundbp[29]);

void performCrankshaftMoves(double stepsize_a, double beta,
    DNAState state[N], DNAState trial[N],
    double BPSEnergies[N-1], double BPSEnergiesTrial[N-1],
    double BPSEnergies2[N-1], double BPSEnergiesTrial2[N-1],
    double p, double *e_A, double *e_B, double *H_e,
    long int *tries_crank, long int *accept_crank, gsl_rng *r, int boundbp[29]);



void perform_cluster_trans_moves(double stepsize_t, double beta, double p,
                                    DNAState state[N], DNAState trial[N],
                                    double BPSEnergies[N-1], double BPSEnergiesTrial[N-1],
                                    double BPSEnergies2[N-1], double BPSEnergiesTrial2[N-1],
                                    double *e_A, double *e_B, double *H_e,
                                    long int *tries_cluster, long int *accept_cluster, gsl_rng *r, int boundbp[29]);


void printReplicas(Replica **replicas, int num_replicas);

void mc_steps(Replica *rep);
void mc_steps_shuffle(Replica *rep);
void shuffle_array(int *array, int n, gsl_rng *r);

// Energy functions
double BPSEnergy(DNAState BP1, DNAState BP2, Base B1, Base B2);
double Energy(DNAState state[N], Base sequence[N]);


// Auxiliary functions
// Function to compute bound sites
void compute_bound_sites(int left, int right, int boundbp[29]);

//// Related to bound base pairs in the nucleosome
int isBound(int base, int boundbp[29]);
int BaseToBond(int base, int boundbp[29]);

//// Related to midplanes
void midplane(DNAState BP1, DNAState BP2, DNAState* out);
void midplanestate(DNAState state[N], DNAState midstate[N-1]);

//// Matrix manipulations
void mmult2(double A[3][3], double B[3][3], double C[3][3]);
void LmmultT2(double A[3][3], double B[3][3], double C[3][3]);
void RmmultT2(double A[3][3], double B[3][3], double C[3][3]);
void OrthogonalizeR(double R[3][3]);
void mat_vec_mult(double A[3][3], double x[3], double out[3]);
double dot(const double a[3], const double b[3]);
double vec_norm(double v[3]);

//// Numerical calculations
void randUnitVectorB(double n[3], gsl_rng *r);
int LinRandBP(int bp, double a);
void RfromAA(double angle, double axis[3], double R[3][3]);

//// Outputting and copying data structures
void PrintBPState(FILE* output, DNAState state);
void PrintDNAState(FILE* output, DNAState state[N]);
void CopyDNAState(DNAState state1[N], DNAState state2[N], int start, int end);
void ReadDNAState(FILE* file, DNAState state[N]);
void CopyBPSEnergies(double Esrc[N-1], double Etrg[N-1], int start, int end);
void PrintSequence(FILE* file, Base* sequence);
void saveStateToPDB(const char *filename, DNAState state[N]);
void saveStateToXYZ(const char *filename, DNAState state[N]);
void combineXYZFiles();	 
void cranklength(int m , int *id_S, int *id_E, gsl_rng *r, int boundbp[29]);
void PrintBPState_terminal(DNAState state);


void segment_rotation(DNAState conf[N], double Rlab[3][3], int idA, int hingedist);
void segment_rotation2(DNAState conf[N], double Rlab[3][3], int id_S, int hingedist);
void check_orthogonal(double R[3][3]);
double realtive_angle(DNAState BP1, DNAState BP2);
Base CharToBase(char c);


/////Parallel Tempering

void TEMPcopyState(DNAState *dest, const DNAState *src, int L);
Replica* createReplica();
void TEMPcopyEnergies(double *dest, const double *src, int L);