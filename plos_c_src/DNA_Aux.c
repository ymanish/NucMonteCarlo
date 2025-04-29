#include "dnaMC2.h"

// Print Base Pair State Configuration (PrintBPState)
//
// Prints a single bp DNAState to the specified output 
// filestream (can also be stdout/stderr). The format is x, y, 
// z, R00, R01, R02, etc., tab-separate, no trailing tab. The 
// printing commands are separated out purely for readability.
//
// Return value: 
// none
//
// Arguments: 
// output - output filestream
// state  - single DNAState instance representing one base pair

void check_orthogonal(double R[3][3]) {

    int i, j, k;
    double tol = 1e-8;

    double RtR[3][3] = { {0.0, 0.0, 0.0},
                          {0.0, 0.0, 0.0},
                          {0.0, 0.0, 0.0} };

    // Compute R^T * R.
    // Note: The (i,j)th entry of R^T * R is sum_{k} R[k][i] * R[k][j]
    for (i = 0; i < 3; i++){
        for (j = 0; j < 3; j++){
            for (k = 0; k < 3; k++){
                RtR[i][j] += R[k][i] * R[k][j];
            }
        }
    }
    
    // Compute the Frobenius norm of (RtR - I)
    double error = 0.0;
    for (i = 0; i < 3; i++){
        for (j = 0; j < 3; j++){
            double diff = RtR[i][j] - ((i == j) ? 1.0 : 0.0);
            error += diff * diff;
        }
    }
    error = sqrt(error);
    
    // If the error is less than tol, we consider R to be orthogonal.
    int value_flag = (error < tol) ? 1 : 0;

    if (value_flag == 1)
        printf("The matrix R is orthogonal.\n");
    else
        printf("The matrix R is NOT orthogonal.\n");


 }







void combineXYZFiles() {
    // Command to combine all .xyz files into a single trajectory file
    const char *command = "cat ./Trajectories/Frame_*.xyz > NucleosomeTrajectory.xyz";

    // Execute the command
    int result = system(command);

    // Check if the command executed successfully
    if (result == 0) {
        printf("Successfully combined .xyz files into NucleosomeTrajectory.xyz\n");
    } else {
        fprintf(stderr, "Error: Failed to combine .xyz files\n");
    }

}

void PrintBPState(FILE* output, DNAState state) {
    // Print position vector components
    fprintf(output, "%e\t%e\t%e\t", state.pos[0], state.pos[1], state.pos[2]);
    // Print orientation matrix components
    fprintf(output, "%e\t%e\t%e\t", state.R[0][0], state.R[0][1], state.R[0][2]);
    fprintf(output, "%e\t%e\t%e\t", state.R[1][0], state.R[1][1], state.R[1][2]);
    fprintf(output, "%e\t%e\t%e", state.R[2][0], state.R[2][1], state.R[2][2]);
}


void PrintBPState_terminal(DNAState state) {
    printf("Position: %f, %f, %f\n", state.pos[0], state.pos[1], state.pos[2]);

    printf("Orientation:\n");
    printf("%f, %f, %f\n", state.R[0][0], state.R[0][1], state.R[0][2]);
    printf("%f, %f, %f\n", state.R[1][0], state.R[1][1], state.R[1][2]);
    printf("%f, %f, %f\n", state.R[2][0], state.R[2][1], state.R[2][2]);
}
// Print DNA Chain Configuration (PrintDNAState)
//
// Calls PrintBPState for every base pair in a DNAState
// array, with a tab character in between.
//
// Return value:
// none
//
// Arguments:
// output - output filestream
// state  - array of DNAStates representing a DNA molecule
void PrintDNAState(FILE* output, DNAState state[N]) {
    int n;
    // For every base pair...
    for ( n = 0; n < N; n++ ) {
        // ...print its state.
        PrintBPState(output, state[n]);
        // Separate by tab characters, but don't add
        // a trailing tab.
        if ( n == N-1 ) { }
        else { fprintf(output, "\t"); }
    }
    fprintf(output, "\n");
}

// Copy DNA Chain Configuration (CopyDNAState)
//
// Copies (a subset of) the base pair configurations in
// a DNAState array to another array. Allows copying only
// a part of the configuration. This allows the function
// to be more efficient when copying trial states to current 
// states in the MC simulation when the MC move only affects
// a small part of the chain.
//
// Return value:
// none
//
// Arguments:
// state1 - source array
// state2 - target array
// start  - base pair index from which to start copying
// stop   - base pair index at which to stop copying
void CopyDNAState(DNAState state1[N], DNAState state2[N], int start, int end) {
    int n;
    // These checks can be optimized away if you're sure you're
    // calling the function correctly.
    
    // In case the start and end arguments are accidentally reversed.
    if ( end < start ) { int x = end; end = start; start = x; }
    
    // In case the start and end are not in the valid range.
    if ( start < 0 ) { start = 0; }
    if ( end > N-1 ) { end = N-1; }
    
    // Copy.
    for ( n = start; n <= end; n++ ) {
        state2[n] = state1[n];
    }
}

// Copy Base Pair Step Elastic Energies (CopyBPSEnergies)
//
// Copies (a subset of) the base pair step energies in
// an array to another array. Its use is similar to that of
// CopyDNAState, but it operates on the BPSEnergies(Trial)
// bookkeeping arrays.
//
// Return value:
// none
//
// Arguments:
// Esrc  - source array
// Etrg  - target array
// start - base pair index from which to start copying
// stop  - base pair index at which to stop copying
void CopyBPSEnergies(double Esrc[N-1], double Etrg[N-1], int start, int end) {
    int n;
    // These checks can be optimized away if you're sure you're
    // calling the function correctly.
    
    // In case the start and end arguments are accidentally reversed.
    if ( end < start ) { int x = end; end = start; start = x; }
    
    // In case the start and end are not in the valid range.
    if ( start < 0 ) { start = 0; }
    if ( end > N-2 ) { end = N-2; }
    
    // Copy.
    for ( n = start; n <= end; n++ ) {
        Etrg[n] = Esrc[n];
    }
}

// Read DNAState From File (ReadDNAState)
//
// Reads in a DNAState from the given file. It expects a 
// specific format: Nx12 tab-separated floating point numbers
// representing in order, x, y, z, R00, R01, R02, etc. for
// N base pairs consecutively. (Note: this is the format in
// which PrintDNAState outputs a DNAState array.) N must be
// equal to the system size set in the header file.
//
// Return value:
// none
//
// Arguments:
// file  - input filestream
// state - DNAState array into which to read the data
void ReadDNAState(FILE* file, DNAState state[N]) {
    int n = 0;
    // Using fscanf to read things in. The syntax is a bit
    // ugly, but that's how it works.
    while ( fscanf(file, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf", &state[n].pos[0], &state[n].pos[1], &state[n].pos[2], &state[n].R[0][0], &state[n].R[0][1], &state[n].R[0][2], &state[n].R[1][0], &state[n].R[1][1], &state[n].R[1][2], &state[n].R[2][0], &state[n].R[2][1], &state[n].R[2][2]) != EOF ) {
        n++;
    }
    // Check that we read in the right amount of data.
    if ( n != N ) { fprintf(stderr, "Error reading data from file: length mismatch, n = %d.\n", n); exit(1); }
}

// Print DNA Sequence (PrintSequence)
//
// Prints a Base array to the specified filestream.
//
// Return value:
// none
//
// Arguments:
// file     - output filestream
// sequence - sequence array to be output
void PrintSequence(FILE* file, Base* sequence) {
    int i;
    char *SeqString;
    // Allocate memory for the string version of the sequence
	SeqString = malloc((N+1)*sizeof(char));
	SeqString[0] = '\0';
	// Convert each character in the sequence to a string
	// and save in the allocated memory
    for ( i = 0; i < N; i++ ) { strcat(SeqString, BaseString[sequence[i]]); }
    // Print the resulting string.
    fprintf(file, "%s\n", SeqString);
}


// Simplified wrappers for some CBLAS matrix multiplication
// functions.

// Basic Matrix Multiplication:
// Multiplies A with B and stores result in C.
void mmult2(double A[3][3], double B[3][3], double C[3][3])  {
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 3, 3, 1.0, *A, 3, *B, 3, 0.0, *C, 3);
}

// Matrix Multiplication with left matrix transposed:
// Multiplies A^T with B and stores result in C.
void LmmultT2(double A[3][3], double B[3][3], double C[3][3])  {
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, 3, 3, 3, 1.0, *A, 3, *B, 3, 0.0, *C, 3);
}

// Matrix Multiplication with right matrix transposed:
// Multiplies A with B^T and stores result in C.
void RmmultT2(double A[3][3], double B[3][3], double C[3][3])  {
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, 3, 3, 3, 1.0, *A, 3, *B, 3, 0.0, *C, 3);
}


void mat_vec_mult(double A[3][3], double x[3], double out[3]) {
    for (int i = 0; i < 3; i++) {
        out[i] = 0.0;
        for (int j = 0; j < 3; j++) {
            out[i] += A[i][j] * x[j];
        }
    }
}



// Helper: Compute dot product of two 3-vectors.
double dot(const double a[3], const double b[3]) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

// Helper: Normalize a 3-vector in place.
double vec_norm(double v[3]) {
    double norm = sqrt(dot(v, v));
    // if (norm > 0)
    //     for (int i = 0; i < 3; i++) v[i] /= norm;
    return norm;
}

double realtive_angle(DNAState BP1, DNAState BP2) {

    double x[3], R[3][3];
    double rot_ang;
    for (int i = 0; i < 3; i++ ) { ;
        x[i] = 0.0;

    for ( int j = 0; j < 3; j++ ) {
        x[i] += BP1.R[j][i]*(BP2.pos[j] - BP1.pos[j]); 
    }
    }

    LmmultT2(BP1.R, BP2.R, R); 

    rot_ang = acos(0.5*(R[0][0] + R[1][1] + R[2][2] - 1.0));
    return rot_ang;
}

// Calculate BPS Midplane (midplane)
//
// Calculates the midplane position and orientation
// between two base pairs.
//
// Return value:
// none
//
// Arguments:
// BP1 - first base pair DNAState
// BP2 - second base pair DNAState
// out - DNAState object into which to write the result
void midplane(DNAState BP1, DNAState BP2, DNAState *out) {
    double R[3][3], Rm[3][3];
    double a[3];
    double rot_ang, frac, t;
    int i;
    
    // Get the rotation matrix between BP1 and BP2
    RmmultT2(BP2.R, BP1.R, R);
    
    // Get the rotation angle
    rot_ang = acos(0.5*(R[0][0] + R[1][1] + R[2][2] - 1.0));
    
    // Get the rotation axis
    frac = 0.5*rot_ang / sin(rot_ang);
    a[0] = frac*(R[2][1] - R[1][2]);
    a[1] = frac*(R[0][2] - R[2][0]);
    a[2] = frac*(R[1][0] - R[0][1]);
	t = 1.0/sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
	for ( i = 0; i < 3; i++ ) {
	    a[i] *= t;
	}
	
	// Generate rotation matrix from the calculated
	// rotation axis, and half the rotation angle.
	// This rotation matrix rotates from BP1 to the
	// midframe.
	RfromAA(0.5*rot_ang, a, Rm);
	
	// Get the midframe orientation
	mmult2(Rm, BP1.R, (*out).R);
	
	// Get the midframe position
	for ( i = 0; i < 3; i++ ) {
	    (*out).pos[i] = 0.5*(BP1.pos[i] + BP2.pos[i]);
	}
}

// Calculate Midplane State (midplanestate)
// 
// Calculates the midplanes between all successive base pairs
// in a DNA chain state.
//
// Return value:
// none
//
// Arguments:
// 
// state    - DNA molecule state
// midstate - DNAState array into which to write the results
void midplanestate(DNAState state[N], DNAState midstate[N-1]) {
    int n;
    for ( n = 0; n < N-1; n++ ) {
        midplane(state[n], state[n+1], &midstate[n]);
    }
}


// Pick random base pair from lineair distribution around bp
// i.e. bp has the highest probability and away from bp the
// probability scales as 1 - a*|n - bp|/N.
// (This function not in use in the Nucleosome code, but can
// come in handy for crankshaft moves and such. Leaving it in
// in case.)
//
// Return value:
// m  - selected base pair index
//
// Arguments:
//
// bp - base pair index around which to center the distribution
// a  - slope of the linear decay
// int LinRandBP(int bp, double a) {
//     double random, func;
//     int m, accept = 0;
//     while ( accept == 0 ) {
//         // Random float between 0 and 1
//         random = gsl_rng_uniform(r);
        
//         // Random integer between 0 and N-1
//         m = gsl_rng_uniform_int(r, N-1);
        
//         // Linear decay function
//         func = 1.0 - a*fabs((double) (m-bp))/N;
        
//         // Accept the trial with probability
//         // given by func
//         if ( random < func ) { accept = 1; }
//     }
//     return m;
// }


// Generate a random unit vector, such that the
// probability density on the surface of the
// unit sphere is uniform. See 
// mathworld.wolfram.com/SpherePointPicking.html
//
// Return value:
// None
//
// Arguments:
// n - 3-array into which to store the result
void randUnitVectorB(double n[3], gsl_rng *r) {

    double a, x, b;
   
    // Random azimuth
    a = gsl_rng_uniform(r)*2*PI;
    
    // Random float between 0 and 1
    x = gsl_rng_uniform(r);
    
    // Random elevation
    b = acos(1-2*x);
    
    // Generate unit vector
    n[0] = sin(b)*cos(a);
    n[1] = sin(b)*sin(a);
    n[2] = cos(b);
}

// Rotation Matrix from Axis-Angle (RfromAA)
//
// Generates a rotation matrix from a given
// unit vector and angle.
//
// Return value:
// none
//
// Arguments:
// angle - rotation angle
// axis  - rotation axis
// R     - 3x3 array into which to store the result
void RfromAA(double angle, double axis[3], double R[3][3]) {
    double c = cos(angle);
    double s = sin(angle);
    double c1 = 1.0-c;
    double s2 = 2.0*s;
    
    R[0][0] = c + axis[0]*axis[0]*c1;
    
    R[1][0] = axis[0]*axis[1]*c1 + axis[2]*s;
    R[2][0] = axis[0]*axis[2]*c1 - axis[1]*s;
    R[0][1] = R[1][0] - axis[2]*s2;
    
    R[1][1] = c + axis[1]*axis[1]*c1;
    
    R[2][1] = axis[1]*axis[2]*c1 + axis[0]*s;
    R[0][2] = R[2][0] + axis[1]*s2;
    R[1][2] = R[2][1] - axis[0]*s2;
    
    R[2][2] = c + axis[2]*axis[2]*c1;
}

// Convert characters to Base variables (CharToBase)
//
// Given a character (one-letter string), give the
// corresponding Base variable. Character can be
// upper or lower case. Any character that is not
// A/a, T/t, C/c or G/g is converted into X (i.e.
// homogeneous DNA).
Base CharToBase(char c) {
	Base B;
	if ( c == 'A' || c == 'a' ) { B = A; }
	else if ( c == 'T' || c == 't' ) { B = T; }
	else if ( c == 'C' || c == 'c' ) { B = C; }
	else if ( c == 'G' || c == 'g' ) { B = G; }
	else { B = X; }
	return B;
}


///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////


void compute_bound_sites(int left, int right, int boundbp[29]) {
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




int isBound(int base, int boundbp[29]) {    
    int base2 = base + 1;
    int n;
    
    boundbp[28] = base;
    
    for ( n = 0; boundbp[n] != base; n++ );
    
    if ( n == 28 ) {
        
        boundbp[28] = base2;
        for ( n = 0; boundbp[n] != base2; n++ );
        
        
        if ( n == 28 ) { return 0; }
        
        else { return 1; }
    }

    return -1;
    
}


int BaseToBond(int base, int boundbp[29]) {
    int n;
    for ( n = 0; n < 28; n++ ) {
        if ( boundbp[n] == base ) {
            return n;
        }
    }
    return -1;
}


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


void cranklength(int m , int *id_S, int *id_E, gsl_rng *r, int boundbp[29]) {


    int mdfrom = boundbp[m]+1;
    int mdto = boundbp[m+1]-1;
    int hinge_dist = 0;
    *id_S = 0;
    *id_E = 0;


    if (mdto - mdfrom >= 3) {

        if (mdto - mdfrom == 3) {
            *id_S = mdfrom;
            *id_E = mdto;

        }
        else {
            
            *id_S = mdfrom + gsl_rng_uniform_int(r, (mdto - mdfrom)-2);
            
            hinge_dist = 3 + gsl_rng_uniform_int(r, (mdto - *id_S -2));

            *id_E = *id_S + hinge_dist;
            if (*id_E > mdto) {
                *id_E = mdto;
            }

            }
        
    }
}
