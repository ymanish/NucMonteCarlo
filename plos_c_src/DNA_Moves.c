#include "dnaMC2.h"

// Single Base Pair Move (SBPMove)
//
// This function performs a random move (translation and rotation)
// on a single base pair.
//
// Return value:
// dE               - Energy difference between trial state and current state
//
// Arguments:
// base             - (zero-indexed) number of the base pair to be moved
// stepsize_t       - maximum step allowed in translational DOF
// stepsize_a       - maximum step allowed in rotational DOF
// currentState     - current accepted state of the system
// newState         - DNAState structure in which to store the new state after the move
// BPSEnergies      - Bookkeeping array containing the current elastic energies between successive base pairs
// BPSEnergiesTrial - Bookkeeping array in which to store the modified elastic energies between successive base pairs

// double SBPMove(int base, double stepsize_t, double stepsize_a, DNAState currentState[N], DNAState newState[N], double BPSEnergies[N-1], double BPSEnergiesTrial[N-1]) {
//     // Declaring variables
//     int i;
//     double dE = 0.0;
//     double step;
    
//     // For each of the translational degrees of freedom...
//     for ( i = 0; i < 3; i++ ) {
//         // ...generate a random value between -stepsize_t and +stepsize_t...
//         step = (2.0*gsl_rng_uniform(r)-1.0)*stepsize_t;
//         // ...and add this to the current position of the base pair, storing the result in newState.
//         newState[base].pos[i] = currentState[base].pos[i] + step;
//     }    
    
//     // Declaring variables for random rotation:
    
//     // Random rotation axis
//     double n[3];
//     randUnitVectorB(n);

                   
//     // Random rotation angle
//     double theta = gsl_rng_uniform(r)*stepsize_a;         
    
//     // Obtain rotation matrix            
//     double R[3][3];
//     RfromAA(theta, n, R);           
    
//     // Rotate base pair, storing the result in newState
//     mmult2(R, currentState[base].R, newState[base].R);      
    
//     // Calculate the energy difference on either side of the base pair
    
//     // Unless the very first base pair was moved...
//     if ( base != 0 )    {     
//         // ...calculate the new energy between the moved base pair and the previous one...
//         BPSEnergiesTrial[base-1] = BPSEnergy(newState[base-1], newState[base], Sequence[base-1], Sequence[base]);
//         // ...and store the difference.
//         dE += BPSEnergiesTrial[base-1] - BPSEnergies[base-1];    
//     }
    
//     // Unless the very last base pair was moved...
//     if ( base != N-1 )  {     
//         // ...do the same for the moved base pair and the next one.
//         BPSEnergiesTrial[base] = BPSEnergy(newState[base], newState[base+1], Sequence[base], Sequence[base+1]);    
//         dE += BPSEnergiesTrial[base] - BPSEnergies[base];     
//     }
    
//     // Return the difference in energy between the trial state and the current state.
//     return dE;
// } 

// // Single Base Pair Move, with a Fixed MidPlane (SBPMoveFMP)
// //
// // This function performs a random move (translation and rotation)
// // on a single base pair, but with the constraint that the base pair
// // has a complementary base pair that needs to be moved in the opposite
// // direction in order to preserve the midplane between the two.
// //
// // Return value:
// // dE               - Energy difference between trial state and current state
// //
// // Arguments:
// // base             - (zero-indexed) number of the base pair to be moved
// // complement       - (zero-indexed) number of the complementary base pair
// // stepsize_t       - maximum step allowed in translational DOF
// // stepsize_a       - maximum step allowed in rotational DOF
// // currentState     - current accepted state of the system
// // newState         - DNAState structure in which to store the new state after the move
// // BPSEnergies      - Bookkeeping array containing the current elastic energies between successive base pairs
// // BPSEnergiesTrial - Bookkeeping array in which to store the modified elastic energies between successive base pairs
// // midplanes        - Array of the orientations and locations of the 28 midplanes that are fixed in the nucleosome

// double SBPMoveFMP(int base, int complement, double stepsize_t, double stepsize_a, 
// DNAState currentState[N], DNAState newState[N], double BPSEnergies[N-1], double BPSEnergiesTrial[N-1],
//  DNAState midplanes[28]) {
//     // Declaring variables
//     int i;
//     double dE = 0.0;
//     double step;
    
//     // Translational move identical to that in SBPMove(...), but the complement is
//     // moved by the same amount in the opposite direction.
//     for ( i = 0; i < 3; i++ ) {
//         step = (2.0*gsl_rng_uniform(r)-1.0)*stepsize_t;
//         newState[base].pos[i] = currentState[base].pos[i] + step;
//         newState[complement].pos[i] = currentState[complement].pos[i] - step;
//     }

//     // The implementation of the random rotation is somewhat non-trivial. One could just generate
//     // a random rotation matrix and apply it to one base pair, while applying its inverse to the
//     // complementary base pair. However, this approach could lead to a numerically unstable midplane.
//     // To ensure that the midplane is absolutely fixed, we instead do the following:
//     //  * Generate random rotation R
//     //  * Calculate the rotation matrix Rmid between the orientation of the base pair and the midplane
//     //  * Calculate the combined rotation matrix R.Rmid, i.e. the rotation between base pair and midplane plus a random rotation
//     //  * Apply this to the midplane orientation to generate a randomly modified orientation for the base pair
//     //  * Apply the inverse rotation to the complementary base pair to generate the appropriate orientation for
//     //    the complementary base pair, regardless of its original orientation.
//     // This scheme will rigidly enforce that the simulation always respects the original midframes.
    
//     // Random rotation matrix
//     double n[3];
//     double R[3][3];
//     double theta = (2.0*gsl_rng_uniform(r)-1.0)*stepsize_a; 
//     randUnitVectorB(n);
//     RfromAA(theta, n, R);

//     // Here we find out which bond we are dealing with, in order to get the right midplane.
//     // Note that we assume that the function has been passed a correct set of base pairs.
//     // ToDo: Error Checking on this point.
//     int bond = BaseToBond(base);
    
//     // If the base pair is not in the list of bound base pairs, we assume the other one must be.
//     if ( bond == -1 ) { bond = BaseToBond(complement); }
        
//     // We obtain the rotation matrix Rmid between the base pair and the midplane...
//     double Rmid[3][3], Rmove[3][3];
//     RmmultT2(currentState[base].R, midplanes[bond].R, Rmid);
    
//     // ...reorthogonalize it (because all the matrix multiplications
//     // tend to be numerical unstable)...
//     OrthogonalizeR(Rmid);

//     // ...and combine it with our random matrix to obtain our move relative to the midframe
//     mmult2(R, Rmid, Rmove);
    
//     // Rotate one base pair by Rmove... 
//     cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 3, 3, 1.0, *Rmove, 3, *midplanes[bond].R, 3, 0.0, *newState[base].R, 3);
//     // ...and the other by the inverse of Rmove.
//     cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, 3, 3, 3, 1.0, *Rmove, 3, *midplanes[bond].R, 3, 0.0, *newState[complement].R, 3);
    
//     // Calculate the new energies and the energy difference
    
//     BPSEnergiesTrial[base-1] = BPSEnergy(newState[base-1], newState[base], Sequence[base-1], Sequence[base]);    
//     dE += BPSEnergiesTrial[base-1] - BPSEnergies[base-1];

//     BPSEnergiesTrial[complement] = BPSEnergy(newState[complement], newState[complement+1], Sequence[complement], Sequence[complement+1]);
//     dE += BPSEnergiesTrial[complement] - BPSEnergies[complement]; 
                 
//     BPSEnergiesTrial[base] = BPSEnergy(newState[base], newState[complement], Sequence[base], Sequence[complement]);       
//     dE += BPSEnergiesTrial[base] - BPSEnergies[base];
    
//     return dE;
// } 

// // Mutation Move
// //
// // This function performs a random mutation on a single base pair.
// //
// // Return value:
// // dE               - Energy difference between trial state and current state
// //
// // Arguments:
// // state            - current state of the simulation
// // seq              - current sequence of the simulation
// // a                - (zero-indexed) number of the base to mutate
// // b                - Base variable in which to store the trial mutation
// // BPSEnergies      - Bookkeeping array containing the current elastic energies between successive base pairs
// // BPSEnergiesTrial - Bookkeeping array in which to store the modified elastic energies between successive base pairs

// double Mutation(DNAState state[N], Base* seq, int a, Base* b, double BPSEnergies[N-1], double BPSEnergiesTrial[N-1]) {
//     // Declare variables
// 	double dE = 0.0;
// 	double s;
	
// 	// Randomly pick a new base identity. The mutation must constitute an actual alteration, hence the while loop.
// 	*b = seq[a];
//     while ( *b == seq[a] ) {
//         s = gsl_rng_uniform(r);
//         if ( s < 0.25 ) { *b = C; }
//         else if ( s < 0.50 ) { *b = T; }
//         else if ( s < 0.75 ) { *b = A; }
//         else { *b = G; }
//     }
	
// 	// Calculate the energy change
//     if ( a > 0 ) {      BPSEnergiesTrial[a-1] = BPSEnergy(state[a-1], state[a], seq[a-1], b[0]);
//                         dE += BPSEnergiesTrial[a-1] - BPSEnergies[a-1];  }
//     if ( a < N-1 ) {    BPSEnergiesTrial[a] = BPSEnergy(state[a], state[a+1], b[0], seq[a+1]);
//                         dE += BPSEnergiesTrial[a] - BPSEnergies[a];  }
        
//     return dE;
// }



// void segment_rotation(DNAState conf[N], double Rlab[3][3], int idA, int hingedist) {
//     int id, idm1;
//     double diff[3], temp_out[3], new_R[3][3], new_pos[3];
//     double rotated_vecs[hingedist - 2][3]; 
    
    

//     id = idA;
//     for (int i = 0; i < hingedist - 2; i++) {
//         idm1 = id;
//         id = (id + 1);
//         for (int j=0; j<3; j++) {diff[j]= conf[id].pos[j] - conf[idm1].pos[j];} // diff = conf[id].pos - conf[idm1].pos
//         // rotated_vec = Rlab * diff
//         cblas_dgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1.0, *Rlab, 3, diff, 1, 0.0, temp_out, 1);

//         // mat_vec_mult(Rlab, diff, temp_out); 
//         for (int j=0; j<3; j++) {rotated_vecs[i][j] = temp_out[j];} 

//     }
    
//     // Then, rebuild positions and update orientations for intermediate base pairs.
//     id = idA;
//     for (int i = 0; i < hingedist - 2; i++) {
//         idm1 = id;
//         id = (id + 1);

//         for (int k = 0; k < 3; k++) {
//         conf[id].pos[k] = conf[idm1].pos[k]+rotated_vecs[i][k];
//         }
//         // Update orientation: new orientation = Rlab * (old orientation)
        
//         mmult2(Rlab, conf[id].R, new_R);
        
//         for (int a = 0; a < 3; a++)
//             for (int b = 0; b < 3; b++)
//                 conf[id].R[a][b] = new_R[a][b];
//     }
// }






// ////////////////////////////////////

// double CrankShaftMove(int idA, int idB, double stepsize_a, double b, DNAState currentState[N], DNAState newState[N], double BPSEnergies[N-1], double BPSEnergiesTrial[N-1]) {
//     // Declaring variables
//     int i,j;
//     double dE = 0.0;
//     double step;
//     int hinge_dist = idB - idA;
//     int idBn = idB -1;
//     int idAn = idA -1;
//     double e, d, w, s;

//     double disp[3] = {0.0, 0.0, 0.0};
//     double rotated_disp[3] = {0.0, 0.0, 0.0};

//     // Declaring variables for random rotation:
    
//     // Random rotation axis
//     double n[3];
//     // randUnitVectorB(n);
//     // Random rotation angle
//     double theta = (2*gsl_rng_uniform(r)-1)*stepsize_a;         

//     double diff[3];

//     for (i = 0; i < 3; i++) {
//         diff[i] = currentState[idB].pos[i] - currentState[idA].pos[i];
//     }
//     // Normalize diff and scale by theta to get the rotation vector Theta.
//     double norm_diff = sqrt(diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2]);

//     for (i = 0; i < 3; i++) {
//         n[i] = diff[i] / norm_diff;
//     }
//     // printf("Normal Vector........... %f %f %f\n", n[0], n[1], n[2]);
//     // Convert Theta (axis-angle) into rotation matrix Rlab.
//     double Rlab[3][3];
//     RfromAA(theta, n, Rlab);

//     // check_orthogonal(Rlab);

//     // ----- Compute Trial Configurations for Boundaries -----
//     // Lower boundary: copy configuration at idA and update its rotation.
//     DNAState TA_rot = currentState[idA];

//     double tempR[3][3];

//     mmult2(Rlab, currentState[idA].R, tempR);
//     for (i = 0; i < 3; i++){for (j = 0; j < 3; j++){TA_rot.R[i][j] = tempR[i][j];}}
        
//     // (Position remains unchanged for TA_rot.)

//     // Upper boundary: update configuration for idBn.
//     DNAState TBn_rot = currentState[idBn];

//     mmult2(Rlab, currentState[idBn].R, tempR);

//     for (i = 0; i < 3; i++){for (j = 0; j < 3; j++){ TBn_rot.R[i][j] = tempR[i][j];}}
        

//     // Update TBn_rot's position:
//     // new position = position(idB) - Rlab * ( position(idB) - position(idBn) )
//     double diff2[3];
//     for (i = 0; i < 3; i++) {
//         diff2[i] = currentState[idB].pos[i] - currentState[idBn].pos[i];
//     }

//     double rotated_diff2[3] = {0.0, 0.0, 0.0};

//     cblas_dgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1.0, *Rlab, 3, diff2, 1, 0.0, rotated_diff2, 1);

//     for (i = 0; i < 3; i++) {
//         TBn_rot.pos[i] = currentState[idB].pos[i] - rotated_diff2[i];
//     }

//     newState[idA] = TA_rot;
//     newState[idBn] = TBn_rot;

//     // mmult2(Rlab, currentState[idB].R, tempR);

//     // for (i = 0; i < 3; i++){for (j = 0; j < 3; j++){currentState[idB].R[i][j] = tempR[i][j];}}

//     // printf("Before rotation Norm of TB - TBn : %f\n", sqrt(diff2[0]*diff2[0] + diff2[1]*diff2[1] + diff2[2]*diff2[2]));

//     // printf("After rotation Norm of TB - TBn : %f\n", sqrt(rotated_diff2[0]*rotated_diff2[0] + rotated_diff2[1]*rotated_diff2[1] + rotated_diff2[2]*rotated_diff2[2]));



//     BPSEnergiesTrial[idAn] = BPSEnergy(currentState[idAn], TA_rot, Sequence[idAn], Sequence[idA]);

//     dE += BPSEnergiesTrial[idAn] - BPSEnergies[idAn];


//     BPSEnergiesTrial[idBn] = BPSEnergy(TBn_rot, currentState[idB], Sequence[idBn], Sequence[idB]);

//     dE += BPSEnergiesTrial[idBn] - BPSEnergies[idBn];


//     if (dE < 0) {
        
//         // currentState[idA] = TA_rot;
//         // currentState[idBn] = TBn_rot;
//         //Accept the move and rotate the rest of the chain.

//     //    printf("Accepted Crankshaft Move: Rotate the rest of the chain\n");

//         //  for (int k = idA; k <= idB-1; k++) {for (i = 0; i < 3; i++) {diff[i] = newState[k+1].pos[i] - newState[k].pos[i];}
//         //     printf("ID %d  norm_diff %f\n", k, vec_norm(diff));
//         //     printf("Relative Angle %f\n", realtive_angle(newState[k], newState[k+1]));

//         // }

//         segment_rotation(newState, Rlab, idA, hinge_dist);

//         // printf("Accepted Crankshaft Move\n");
//         // printf("Checking Orthogonality of newstates\n");

//         // for (int idx = idA; idx <= idB; idx++) {
//         //     check_orthogonal(newState[idx].R);
//         // }


//         // for (int k = idA; k <= idB-1; k++) {for (i = 0; i < 3; i++) {diff[i] = newState[k+1].pos[i] - newState[k].pos[i];}
//         //     printf("ID %d  norm_diff %f\n", k, vec_norm(diff));
//         //     printf("Relative Angle %f\n", realtive_angle(newState[k], newState[k+1]));
//         // }


//         return dE;


//     } else {
//         // Random float between 0 and 1.
//         d = gsl_rng_uniform(r);
        
//         w = exp(-b*dE); 
        
//         // With probability given by w, accept the move anyway.
//         if ( d < w ) {
            

//             // printf("Accepted Crankshaft Move: Rotate the rest of the chain\n");
//             // for (int k = idA; k <= idB-1; k++) {for (i = 0; i < 3; i++) {diff[i] = newState[k+1].pos[i] - newState[k].pos[i];}
//             //     printf("ID %d  norm_diff %f\n", k, vec_norm(diff));
//             //     printf("Relative Angle %f\n", realtive_angle(newState[k], newState[k+1]));
//             // }


//             // for (int idx = idA; idx <= idB; idx++) {
//             //     check_orthogonal(newState[idx].R);
//             //     PrintBPState_terminal(newState[idx]);

//             // }

//             segment_rotation(newState, Rlab, idA, hinge_dist);

//             // printf("Accepted Crankshaft Move\n");
//             // printf("Checking Orthogonality of newstates\n");

//             // for (int idx = idA; idx <= idB; idx++) {
//             //     check_orthogonal(newState[idx].R);
//             //     PrintBPState_terminal(newState[idx]);

//             // }

//             // for (int k = idA; k <= idB-1; k++) {for (i = 0; i < 3; i++) {diff[i] = newState[k+1].pos[i] - newState[k].pos[i];}
//             //     printf("ID %d  norm_diff %f\n", k, vec_norm(diff));
//             //     printf("Relative Angle %f\n", realtive_angle(newState[k], newState[k+1]));
//             // }



//             return dE; 

//         } 
//         else {
//             // Reject the move.
//             // printf("Rejected Crankshaft Move\n");
//             dE = 0.0;
            
//             newState[idA] = currentState[idA];
//             newState[idBn] = currentState[idBn];

//             BPSEnergiesTrial[idAn] = BPSEnergies[idAn];
//             BPSEnergiesTrial[idBn] = BPSEnergies[idBn];

//             return dE;
//         }
//     }
                   
// } 








double SBPMove2(int base, double stepsize_t, double stepsize_a, DNAState currentState[N], DNAState newState[N], double BPSEnergies[N-1], double BPSEnergiesTrial[N-1], double BPSEnergies2[N-1], double BPSEnergiesTrial2[N-1], double p, double *de_a, double *de_b, gsl_rng *r) {
    // Declaring variables
    int i;
    // double dE = 0.0;
    // double de_a = 0.0;
    // double de_b = 0.0;

    double dH = 0.0; 

    *de_a = 0.0;
    *de_b = 0.0;

    double step;

    step = (2.0*gsl_rng_uniform(r)-1.0)*stepsize_t;

    double nt[3];
    randUnitVectorB(nt, r);

    double ntscaled[3];
    for (i = 0; i < 3; i++) {
        ntscaled[i] = step * nt[i];
    }


    // For each of the translational degrees of freedom...
    for ( i = 0; i < 3; i++ ) {
        // ...generate a random value between -stepsize_t and +stepsize_t...
        // step = (2.0*gsl_rng_uniform(r)-1.0)*stepsize_t;
        // ...and add this to the current position of the base pair, storing the result in newState.
        newState[base].pos[i] = currentState[base].pos[i] + ntscaled[i];
    }    
    
    // Declaring variables for random rotation:
    
    // Random rotation axis
    double n[3];
    randUnitVectorB(n, r);

    
    // fprintf(stderr, "n: %f %f %f\n", n[0], n[1], n[2]);
    // exit(0);
                   
    // Random rotation angle
    double theta = gsl_rng_uniform(r)*stepsize_a;         
    
    // Obtain rotation matrix            
    double R[3][3];
    RfromAA(theta, n, R);           
    
    // Rotate base pair, storing the result in newState
    mmult2(R, currentState[base].R, newState[base].R);      
    
    // Calculate the energy difference on either side of the base pair
    
    // Unless the very first base pair was moved...
    if ( base != 0 )    {     
        // ...calculate the new energy between the moved base pair and the previous one...
        BPSEnergiesTrial[base-1] = BPSEnergy(newState[base-1], newState[base], Sequence[base-1], Sequence[base]);
        // ...and store the difference.
        *de_a += BPSEnergiesTrial[base-1] - BPSEnergies[base-1];    

        BPSEnergiesTrial2[base-1] = BPSEnergy(newState[base-1], newState[base], Sequence2[base-1], Sequence2[base]);
        // ...and store the difference.
        *de_b += BPSEnergiesTrial2[base-1] - BPSEnergies2[base-1];
    }
    
    // Unless the very last base pair was moved...
    if ( base != N-1 )  {     
        // ...do the same for the moved base pair and the next one.
        BPSEnergiesTrial[base] = BPSEnergy(newState[base], newState[base+1], Sequence[base], Sequence[base+1]);    
        *de_a += BPSEnergiesTrial[base] - BPSEnergies[base];    

        BPSEnergiesTrial2[base] = BPSEnergy(newState[base], newState[base+1], Sequence2[base], Sequence2[base+1]);
        *de_b += BPSEnergiesTrial2[base] - BPSEnergies2[base]; 
    }
    dH = (1-p)*(*de_a) + p*(*de_b);

    // Return the difference in energy between the trial state and the current state.
    return dH;
} 


double SBPMoveFMP2(int base, int complement, double stepsize_t, double stepsize_a, 
DNAState currentState[N], DNAState newState[N], double BPSEnergies[N-1], double BPSEnergiesTrial[N-1], double BPSEnergies2[N-1], double BPSEnergiesTrial2[N-1],
 DNAState midplanes[28], double p, double *de_a, double *de_b, gsl_rng *r, int boundbp[29]) {
    // Declaring variables
    int i;
    // // double dE = 0.0;
    *de_a = 0.0;
    *de_b = 0.0;
    double dH = 0.0;
    double step;
    
 

    for ( i = 0; i < 3; i++ ) {
        step = (2.0*gsl_rng_uniform(r)-1.0)*stepsize_t;
        newState[base].pos[i] = currentState[base].pos[i] + step;
        newState[complement].pos[i] = currentState[complement].pos[i] - step;
    }

   
    // Random rotation matrix
    double n[3];
    double R[3][3];
    double theta = (2.0*gsl_rng_uniform(r)-1.0)*stepsize_a; 

    randUnitVectorB(n, r);
    RfromAA(theta, n, R);

    int bond = BaseToBond(base, boundbp);
    
    // If the base pair is not in the list of bound base pairs, we assume the other one must be.
    if ( bond == -1 ) { bond = BaseToBond(complement, boundbp); }
        
    // We obtain the rotation matrix Rmid between the base pair and the midplane...
    double Rmid[3][3], Rmove[3][3];
    RmmultT2(currentState[base].R, midplanes[bond].R, Rmid);
    
   
    OrthogonalizeR(Rmid);

    // ...and combine it with our random matrix to obtain our move relative to the midframe
    mmult2(R, Rmid, Rmove);
    
    // Rotate one base pair by Rmove... 
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 3, 3, 1.0, *Rmove, 3, *midplanes[bond].R, 3, 0.0, *newState[base].R, 3);
    // ...and the other by the inverse of Rmove.
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, 3, 3, 3, 1.0, *Rmove, 3, *midplanes[bond].R, 3, 0.0, *newState[complement].R, 3);
    
    // Calculate the new energies and the energy difference
    
    BPSEnergiesTrial[base-1] = BPSEnergy(newState[base-1], newState[base], Sequence[base-1], Sequence[base]);    
    *de_a += BPSEnergiesTrial[base-1] - BPSEnergies[base-1];

    BPSEnergiesTrial[complement] = BPSEnergy(newState[complement], newState[complement+1], Sequence[complement], Sequence[complement+1]);
    *de_a += BPSEnergiesTrial[complement] - BPSEnergies[complement]; 
                 
    BPSEnergiesTrial[base] = BPSEnergy(newState[base], newState[complement], Sequence[base], Sequence[complement]);       
    *de_a += BPSEnergiesTrial[base] - BPSEnergies[base];



    BPSEnergiesTrial2[base-1] = BPSEnergy(newState[base-1], newState[base], Sequence2[base-1], Sequence2[base]);
    *de_b += BPSEnergiesTrial2[base-1] - BPSEnergies2[base-1];

    BPSEnergiesTrial2[complement] = BPSEnergy(newState[complement], newState[complement+1], Sequence2[complement], Sequence2[complement+1]);
    *de_b += BPSEnergiesTrial2[complement] - BPSEnergies2[complement];

    BPSEnergiesTrial2[base] = BPSEnergy(newState[base], newState[complement], Sequence2[base], Sequence2[complement]);
    *de_b += BPSEnergiesTrial2[base] - BPSEnergies2[base];

    // fprintf(stderr, "de_a: %f, de_b: %f\n", *de_a, *de_b);
    dH = (1-p)*(*de_a) + p*(*de_b);

    //  if (isnan(dH)) {
    //         fprintf(stderr, "Inside the Function SBPMoveFMP2 %f \n", step);

    //         printf("stepsize_a = %f\n", stepsize_a);
    //         printf("stepsize_t = %f\n", stepsize_t);
    //         printf("&rep->energy_A = %f\n",*de_a );
    //         printf("&rep->energy_B = %f\n",*de_b);
    //         printf("&rep->H_e = %f\n", dH);
    //         printf("base = %d\n", base);

    //         printf("BPS Energies Trial base -1 %f\n", BPSEnergiesTrial[base-1]);
    //         printf("BPS Energies base -1 %f\n", BPSEnergies[base-1]);
    //         printf("BPS Energies Trial complement %f\n", BPSEnergiesTrial[complement]);
    //         printf("BPS Energies complement %f\n", BPSEnergies[complement]);
    //         printf("BPS Energies Trial base %f\n", BPSEnergiesTrial[base]);
    //         printf("BPS Energies base %f\n", BPSEnergies[base]);


    //         printf("BPS Energies Trial 2 base -1 %f\n", BPSEnergiesTrial2[base-1]);
    //         printf("BPS Energies 2 base -1 %f\n", BPSEnergies2[base-1]);
    //         printf("BPS Energies Trial 2 complement %f\n", BPSEnergiesTrial2[complement]);
    //         printf("BPS Energies 2 complement%f\n", BPSEnergies2[complement]);
    //         printf("BPS Energies Trial 2 base %f\n", BPSEnergiesTrial2[base]);
    //         printf("BPS Energies 2 base %f\n", BPSEnergies2[base]);

    //         printf("Sequence base -1 %d\n", Sequence[base-1]);
    //         printf("Sequence base %d\n", Sequence[base]);
    //         printf("Sequence complement %d\n", Sequence[complement]);

    //         exit(1);
    //          }

    return dH;
} 







void segment_rotation2(DNAState conf[N], double Rlab[3][3], int id_S, int hingedist) {
    int id, idm1;
    double diff[3], temp_out[3], new_R[3][3], new_pos[3];
    double rotated_vecs[hingedist - 2][3]; 
    
    

    id = id_S;
    for (int i = 0; i < hingedist - 2; i++) {
        idm1 = id;
        id = (id + 1);
        for (int j=0; j<3; j++) {diff[j]= conf[id].pos[j] - conf[idm1].pos[j];} // diff = conf[id].pos - conf[idm1].pos
        // rotated_vec = Rlab * diff
        cblas_dgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1.0, *Rlab, 3, diff, 1, 0.0, temp_out, 1);

        // mat_vec_mult(Rlab, diff, temp_out); 
        for (int j=0; j<3; j++) {rotated_vecs[i][j] = temp_out[j];} 

    }
    
    // Then, rebuild positions and update orientations for intermediate base pairs.
    id = id_S;
    for (int i = 0; i < hingedist - 2; i++) {
        idm1 = id;
        id = (id + 1);

        for (int k = 0; k < 3; k++) {
        conf[id].pos[k] = conf[idm1].pos[k]+rotated_vecs[i][k];
        }
        // Update orientation: new orientation = Rlab * (old orientation)
        
        mmult2(Rlab, conf[id].R, new_R);
        
        for (int a = 0; a < 3; a++)
            for (int b = 0; b < 3; b++)
                conf[id].R[a][b] = new_R[a][b];
    }
}






////////////////////////////////////

double CrankShaftMove2(int id_S, int id_E, double stepsize_a, double beta, DNAState currentState[N], DNAState newState[N], 
                        double BPSEnergies[N-1], double BPSEnergiesTrial[N-1],
                         double BPSEnergies2[N-1], double BPSEnergiesTrial2[N-1], 
                         double p, double *de_a, double *de_b, gsl_rng *r) {
    // Declaring variables
    int i,j;
    double dH = 0.0;
    *de_a = 0.0;
    *de_b = 0.0;
    double step;
    int hinge_dist = id_E - id_S;
    int id_En = id_E -1;
    int id_Sn = id_S -1;
    double e, d, w, s;

    double disp[3] = {0.0, 0.0, 0.0};
    double rotated_disp[3] = {0.0, 0.0, 0.0};
    double rotated_diff2[3] = {0.0, 0.0, 0.0};

    double n[3], diff[3], Rlab[3][3], tempR[3][3], diff2[3];

    double theta = (2*gsl_rng_uniform(r)-1)*stepsize_a;         

    for (i = 0; i < 3; i++) {
        diff[i] = currentState[id_E].pos[i] - currentState[id_S].pos[i];
    }
    // Normalize diff and scale by theta to get the rotation vector Theta.
    double norm_diff = sqrt(diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2]);

    for (i = 0; i < 3; i++) {
        n[i] = diff[i] / norm_diff;
    }

    RfromAA(theta, n, Rlab);

    
    DNAState TA_rot = currentState[id_S];

    mmult2(Rlab, currentState[id_S].R, tempR);
    for (i = 0; i < 3; i++){for (j = 0; j < 3; j++){TA_rot.R[i][j] = tempR[i][j];}}
        

    DNAState TBn_rot = currentState[id_En];

    mmult2(Rlab, currentState[id_En].R, tempR);

    for (i = 0; i < 3; i++){for (j = 0; j < 3; j++){ TBn_rot.R[i][j] = tempR[i][j];}}
        

    
    for (i = 0; i < 3; i++) {
        diff2[i] = currentState[id_E].pos[i] - currentState[id_En].pos[i];
    }


    cblas_dgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1.0, *Rlab, 3, diff2, 1, 0.0, rotated_diff2, 1);

    for (i = 0; i < 3; i++) {
        TBn_rot.pos[i] = currentState[id_E].pos[i] - rotated_diff2[i];
    }

    newState[id_S] = TA_rot;
    newState[id_En] = TBn_rot;


    BPSEnergiesTrial[id_Sn] = BPSEnergy(currentState[id_Sn], TA_rot, Sequence[id_Sn], Sequence[id_S]);
    *de_a += BPSEnergiesTrial[id_Sn] - BPSEnergies[id_Sn];
    // printf("Sequence 1 BP_ENERGIES INSIDE THE Cluster Move, idS-1.... %f, %f\n", BPSEnergiesTrial[id_Sn], BPSEnergies[id_Sn]);


    BPSEnergiesTrial[id_En] = BPSEnergy(TBn_rot, currentState[id_E], Sequence[id_En], Sequence[id_E]);
    *de_a += BPSEnergiesTrial[id_En] - BPSEnergies[id_En];
    // printf("Sequence 1 BP_ENERGIES INSIDE THE Cluster Move, idE-1.... %f, %f\n", BPSEnergiesTrial[id_En], BPSEnergies[id_En]);


    BPSEnergiesTrial2[id_Sn] = BPSEnergy(currentState[id_Sn], TA_rot, Sequence2[id_Sn], Sequence2[id_S]);
    *de_b += BPSEnergiesTrial2[id_Sn] - BPSEnergies2[id_Sn];
    // printf("Sequence 2 BP_ENERGIES INSIDE THE Cluster Move, idS-1.... %f, %f\n", BPSEnergiesTrial2[id_Sn], BPSEnergies2[id_Sn]);


    BPSEnergiesTrial2[id_En] = BPSEnergy(TBn_rot, currentState[id_E], Sequence2[id_En], Sequence2[id_E]);
    *de_b += BPSEnergiesTrial2[id_En] - BPSEnergies2[id_En];
    // printf("Sequence 2 BP_ENERGIES INSIDE THE Cluster Move, idE-1.... %f, %f\n", BPSEnergiesTrial2[id_En], BPSEnergies2[id_En]);


    // printf("*de_a: %f, *de_b: %f\n", *de_a, *de_b);
    
    dH = (1-p)*(*de_a) + p*(*de_b);



    if (dH < 0) {
  
        // printf("Accepted Crankshaft Move\n");

        segment_rotation2(newState, Rlab, id_S, hinge_dist);

        return dH;


    } else {
        // Random float between 0 and 1.
        d = gsl_rng_uniform(r);
        
        w = exp(-beta*dH); 
        // printf("dH: %f, boltzmann: %f, random: %f\n", dH, exp(-beta*dH), d);

        // With probability given by w, accept the move anyway.
        if ( d < w ) {
            // printf("Accepted Crankshaft Move\n");
            
            segment_rotation2(newState, Rlab, id_S, hinge_dist);
            // printf("BP Energies: %f, %f, %f, %f\n", BPSEnergiesTrial[id_Sn], BPSEnergiesTrial[id_En], BPSEnergiesTrial2[id_Sn], BPSEnergiesTrial2[id_En]);
            return dH; 

        } 
        else {
            // Reject the move.
            // printf("Rejected Crankshaft Move\n");

            dH = 0.0;
            *de_a = 0.0;
            *de_b = 0.0;
            
            newState[id_S] = currentState[id_S];
            newState[id_En] = currentState[id_En];

            BPSEnergiesTrial[id_Sn] = BPSEnergies[id_Sn];
            BPSEnergiesTrial[id_En] = BPSEnergies[id_En];

            BPSEnergiesTrial2[id_Sn] = BPSEnergies2[id_Sn];
            BPSEnergiesTrial2[id_En] = BPSEnergies2[id_En];


            return dH;
        }
    }
                   
} 



double ClusterTrans2(int id_S, int id_E, double stepsize_t, double beta, DNAState currentState[N], DNAState newState[N], 
                        double BPSEnergies[N-1], double BPSEnergiesTrial[N-1],
                         double BPSEnergies2[N-1], double BPSEnergiesTrial2[N-1], 
                         double p, double *de_a, double *de_b, gsl_rng *r) {

    // printf("Inside the Cluster Move\n");

    int i;
    double dH = 0.0;
    *de_a = 0.0;
    *de_b = 0.0;
    double step;
    // int hinge_dist = id_E - id_S;
    int id_En = id_E -1;
    int id_Sn = id_S -1;
    double e, d, w, s;


    step = (2.0*gsl_rng_uniform(r)-1.0)*stepsize_t;

    // for ( i = 0; i < 3; i++ ) {
    //     newState[base].pos[i] = currentState[base].pos[i] + step;
    // }    
    
    double n[3];
    randUnitVectorB(n, r);

    double nscaled[3];
    for (i = 0; i < 3; i++) {
        nscaled[i] = step * n[i];
    }


    DNAState TA_rot = currentState[id_S];

    DNAState TBn_rot = currentState[id_En];

 
    for (int j = 0; j < 3; j++) { TA_rot.pos[j] = currentState[id_S].pos[j] + nscaled[j]; }
    for (int j = 0; j < 3; j++) { TBn_rot.pos[j] = currentState[id_En].pos[j] + nscaled[j]; }

    // newState[id_S] = TA_rot;
    // newState[id_En] = TBn_rot;


    BPSEnergiesTrial[id_Sn] = BPSEnergy(currentState[id_Sn], TA_rot, Sequence[id_Sn], Sequence[id_S]);
    *de_a += BPSEnergiesTrial[id_Sn] - BPSEnergies[id_Sn];
    // printf("Sequence 1 BP_ENERGIES INSIDE THE Cluster Move, idS-1.... %f, %f\n", BPSEnergiesTrial[id_Sn], BPSEnergies[id_Sn]);


    BPSEnergiesTrial[id_En] = BPSEnergy(TBn_rot, currentState[id_E], Sequence[id_En], Sequence[id_E]);
    *de_a += BPSEnergiesTrial[id_En] - BPSEnergies[id_En];
    // printf("Sequence 1 BP_ENERGIES INSIDE THE Cluster Move, idE-1.... %f, %f\n", BPSEnergiesTrial[id_En], BPSEnergies[id_En]);


    BPSEnergiesTrial2[id_Sn] = BPSEnergy(currentState[id_Sn], TA_rot, Sequence2[id_Sn], Sequence2[id_S]);
    *de_b += BPSEnergiesTrial2[id_Sn] - BPSEnergies2[id_Sn];
    // printf("Sequence 2 BP_ENERGIES INSIDE THE Cluster Move, idS-1.... %f, %f\n", BPSEnergiesTrial2[id_Sn], BPSEnergies2[id_Sn]);


    BPSEnergiesTrial2[id_En] = BPSEnergy(TBn_rot, currentState[id_E], Sequence2[id_En], Sequence2[id_E]);
    *de_b += BPSEnergiesTrial2[id_En] - BPSEnergies2[id_En];
    // printf("Sequence 2 BP_ENERGIES INSIDE THE Cluster Move, idE-1.... %f, %f\n", BPSEnergiesTrial2[id_En], BPSEnergies2[id_En]);


    // printf("*de_a: %f, *de_b: %f\n", *de_a, *de_b);
    
    dH = (1-p)*(*de_a) + p*(*de_b);


    
    if (dH < 0) {
  
        // printf("Accepted Crankshaft Move\n");

        for (int k = id_S; k< id_E; k++) {
            for (int j = 0; j < 3; j++) {
                newState[k].pos[j] = currentState[k].pos[j] + nscaled[j];
            }
        }
        return dH;


    } else {
        // Random float between 0 and 1.
        d = gsl_rng_uniform(r);
        
        w = exp(-beta*dH); 
        // printf("dH: %f, boltzmann: %f, random: %f\n", dH, exp(-beta*dH), d);

        // With probability given by w, accept the move anyway.
        if ( d < w ) {
            // printf("Accepted Crankshaft Move\n");
            
            for (int k = id_S; k< id_E; k++) {
                for (int j = 0; j < 3; j++) {
                newState[k].pos[j] = currentState[k].pos[j] + nscaled[j];
                }
             }           
            // printf("BP Energies: %f, %f, %f, %f\n", BPSEnergiesTrial[id_Sn], BPSEnergiesTrial[id_En], BPSEnergiesTrial2[id_Sn], BPSEnergiesTrial2[id_En]);
            return dH; 

        } 
        else {
            // Reject the move.
            // printf("Rejected Crankshaft Move\n");

            dH = 0.0;
            *de_a = 0.0;
            *de_b = 0.0;
            
            // newState[id_S] = currentState[id_S];
            // newState[id_En] = currentState[id_En];

            BPSEnergiesTrial[id_Sn] = BPSEnergies[id_Sn];
            BPSEnergiesTrial[id_En] = BPSEnergies[id_En];

            BPSEnergiesTrial2[id_Sn] = BPSEnergies2[id_Sn];
            BPSEnergiesTrial2[id_En] = BPSEnergies2[id_En];


            return dH;
        }
    }



}


















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
                        gsl_rng *r, 
                        int boundbp[29]) 
{
    int bb, bplo, bphi;
    double s, d, w, dH;
    double dE_A = 0.0, dE_B = 0.0;

    bb = isBound(n1, boundbp);

    bplo = ( n1 < n1+bb ? n1 : n1+bb ); 
    bphi = ( n1 < n1+bb ? n1+bb : n1 );

    if ( bb ) { 
        (*midbp_tries)++;
        // printf("P value %f\n", p);
        // printf("stepsize_a = %f\n", stepsize_a);
        // printf("stepsize_t = %f\n", stepsize_t);
        // printf("&rep->delta_energy_A = %f\n", dE_A);
        // printf("&rep->delta_energy_B = %f\n",dE_B);
        // printf("&rep->H_e = %f\n", *H_e);
        dH = SBPMoveFMP2(bplo, bphi, stepsize_a, stepsize_t,
                            state, trial, BPSEnergies, BPSEnergiesTrial,
                            BPSEnergies2, BPSEnergiesTrial2, midplanes, p,
                            &dE_A, &dE_B, r, boundbp);



        //  if (isnan(dH)) {
        //     fprintf(stderr, "Error: Energy dH is NaN for replica %f. Inside SBPmOVE\n", beta);
        //     printf("P value %f\n", p);
        //     printf("stepsize_a = %f\n", stepsize_a);
        //     printf("stepsize_t = %f\n", stepsize_t);
        //     printf("&rep->delta_energy_A = %f\n", dE_A);
        //     printf("&rep->delta_energy_B = %f\n",dE_B);
        //     printf("&rep->H_e = %f\n", *H_e);
        //     exit(1);
        //      }
    }
    else {   
        (*sbp_tries)++;  
        dH = SBPMove2(n1, 2.0*stepsize_a, 2.0*stepsize_t,
                        state, trial, BPSEnergies, BPSEnergiesTrial,
                        BPSEnergies2, BPSEnergiesTrial2, p,
                        &dE_A, &dE_B, r);


        //  if (isnan(dH)) {
        //     fprintf(stderr, "Error: Energy dH is NaN in SBPm for replica %f. Inside SBPmOVE\n", beta);

        //     printf("stepsize_a = %f\n", stepsize_a);
        //     printf("stepsize_t = %f\n", stepsize_t);
        //     printf("&rep->energy_A = %f\n", dE_A);
        //     printf("&rep->energy_B = %f\n",dE_B);
        //     printf("&rep->H_e = %f\n", *H_e);
        //     exit(1);
        //      }
    }

    if ( dH < 0 ) { 
        if ( bb ) { 
            (*midbp_accept)++;
        }
        else {
            (*sbp_accept)++;
        }
        CopyDNAState(trial, state, bplo, bphi);
        CopyBPSEnergies(BPSEnergiesTrial, BPSEnergies, bplo-1, bphi);
        CopyBPSEnergies(BPSEnergiesTrial2, BPSEnergies2, bplo-1, bphi);

        *e_A += dE_A;
        *e_B += dE_B;
        *H_e += dH;
    }
    else { 
        d = gsl_rng_uniform(r);
        w = exp(-beta*dH); 
        // printf("dH: %f, boltzmann: %f, random: %f\n", dH, exp(-beta*dH), d);
        if ( d < w ) {
            if ( bb ) { 
                (*midbp_accept)++;
            }
            else {
                (*sbp_accept)++;
            }
            CopyDNAState(trial, state, bplo, bphi);
            CopyBPSEnergies(BPSEnergiesTrial, BPSEnergies, bplo-1, bphi);
            CopyBPSEnergies(BPSEnergiesTrial2, BPSEnergies2, bplo-1, bphi);
            *e_A += dE_A;
            *e_B += dE_B;
            *H_e += dH;
        }
        else {
            CopyDNAState(state, trial, bplo, bphi);
            CopyBPSEnergies(BPSEnergies, BPSEnergiesTrial, bplo-1, bphi);
            CopyBPSEnergies(BPSEnergies2, BPSEnergiesTrial2, bplo-1, bphi);
        }
    }
}







void performCrankshaftMoves(double stepsize_a, double beta,
    DNAState state[N], DNAState trial[N],
    double BPSEnergies[N-1], double BPSEnergiesTrial[N-1],
    double BPSEnergies2[N-1], double BPSEnergiesTrial2[N-1],
    double p, double *e_A, double *e_B, double *H_e,
    long int *tries_crank, long int *accept_crank, gsl_rng *r, int boundbp[29])
{
    int id_S=0;
    int id_E=0;
    double dH=0.0; 
    double dE_A=0.0; 
    double dE_B=0.0;
    
    for (int m = 0; m < 27; m++) {
        cranklength(m, &id_S, &id_E, r, boundbp);
        if (id_S == 0 || id_E == 0) {
            continue;
        } else {
            (*tries_crank)++;

            dH = CrankShaftMove2(id_S, id_E, 15 * stepsize_a, beta,
                state, trial,
                BPSEnergies, BPSEnergiesTrial,
                BPSEnergies2, BPSEnergiesTrial2,
                p, &dE_A, &dE_B, r);

            if (dH != 0.0) {
                (*accept_crank)++;
            }

            CopyDNAState(trial, state, id_S, id_E);
            CopyBPSEnergies(BPSEnergiesTrial, BPSEnergies, id_S - 1, id_E - 1);
            CopyBPSEnergies(BPSEnergiesTrial2, BPSEnergies2, id_S - 1, id_E - 1);

            *e_A += dE_A;
            *e_B += dE_B;
            *H_e += dH;
        }
    }
}




void perform_cluster_trans_moves(double stepsize_t, double beta, double p,
                                    DNAState state[N], DNAState trial[N],
                                    double BPSEnergies[N-1], double BPSEnergiesTrial[N-1],
                                    double BPSEnergies2[N-1], double BPSEnergiesTrial2[N-1],
                                    double *e_A, double *e_B, double *H_e,
                                    long int *tries_cluster, long int *accept_cluster, gsl_rng *r, int boundbp[29]) {
    int id_S=0;
    int id_E=0;
    double dH=0.0; 
    double dE_A=0.0; 
    double dE_B=0.0;

    for (int m = 0; m < 27; m++) {
        cranklength(m, &id_S, &id_E, r, boundbp);

        if (id_S == 0 || id_E == 0)
            continue;
        else {
            (*tries_cluster)++;

            dH = ClusterTrans2(id_S, id_E, 10 * stepsize_t, beta, state, trial, 
                                BPSEnergies, BPSEnergiesTrial, 
                                BPSEnergies2, BPSEnergiesTrial2, p, &dE_A, &dE_B, r);

            if (dH != 0.0)
                (*accept_cluster)++;

            CopyDNAState(trial, state, id_S, id_E);
            CopyBPSEnergies(BPSEnergiesTrial, BPSEnergies, id_S - 1, id_E - 1);
            CopyBPSEnergies(BPSEnergiesTrial2, BPSEnergies2, id_S - 1, id_E - 1);

            *e_A += dE_A;
            *e_B += dE_B;
            *H_e += dH;
        }
    }
}

















////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////   
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




void mc_steps(Replica *rep) {
    

    // // Initialize variables

   
    
    long int sweep_size = 10;
    int C_Max = 1;
    int Local_Max = 10;

    long int sbp_tries = 0; 
    long int midbp_tries = 0;

    long int sbp_accept = 0; 
    long int midbp_accept = 0;

    long int tries_crank = 0;
    long int tries_cluster = 0; 

    long int accept_crank = 0;
    long int accept_cluster = 0; 

    int n1, id_S, id_E, LED_n1, RED_n1;

    double b = rep->beta;
    double stepsize_t = 1.0/(sqrt(b)*10.0);
    double stepsize_a = PI/(sqrt(b)*300.0);
    
    gsl_rng *local_rng = rep->rn;

    

    // compute_bound_sites();
    
    for (int move = 0; move < sweep_size; move++) {

        for (int bp_move = 0; bp_move < Local_Max; bp_move++) {

            n1 = gsl_rng_uniform_int(local_rng, N);

           



            perform_sbp_move(n1, 
                            stepsize_a, 
                            stepsize_t, 
                            rep->beta,
                            rep->state,
                            rep->trial,
                            rep->BPSEnergies,
                            rep->BPSEnergiesTrial,
                            rep->BPSEnergies2,
                            rep->BPSEnergiesTrial2,
                            rep->midplanes,
                            rep->mu,
                            &rep->energy_A,
                            &rep->energy_B,
                            &rep->H_e,
                            &sbp_tries,
                            &midbp_tries,
                            &sbp_accept,
                            &midbp_accept, 
                            local_rng, 
                            rep->boundbp); 

                
            // if (isnan(rep->H_e)) {
            // fprintf(stderr, "Error: Energy H_e is NaN for replica %f. Inside SBPmOVE\n", rep->beta);

            // printf("stepsize_a = %f\n", stepsize_a);
            // printf("stepsize_t = %f\n", stepsize_t);
            // printf("rep->beta = %f\n", rep->beta);
            // printf("&rep->energy_A = %f\n", rep->energy_A);
            // printf("&rep->energy_B = %f\n",rep->energy_B);
            // printf("&rep->H_e = %f\n", rep->H_e);
            // printf("&sbp_tries = %ld\n", sbp_tries);
            // printf("&midbp_tries = %ld\n", midbp_tries);
            // printf("&sbp_accept = %ld\n", sbp_accept);
            // printf("&midbp_accept = %ld\n", midbp_accept);
            // exit(1);
            //  }

                
        }


        LED_n1 = gsl_rng_uniform_int(local_rng, rep->boundbp[0]);
        perform_sbp_move(LED_n1, 
            stepsize_a, 
            stepsize_t, 
            rep->beta,
            rep->state,
            rep->trial,
            rep->BPSEnergies,
            rep->BPSEnergiesTrial,
            rep->BPSEnergies2,
            rep->BPSEnergiesTrial2,
            rep->midplanes,
            rep->mu,
            &rep->energy_A,
            &rep->energy_B,
            &rep->H_e,
            &sbp_tries,
            &midbp_tries,
            &sbp_accept,
            &midbp_accept, local_rng, rep->boundbp);

        //   if (isnan(rep->H_e)) {
        //     fprintf(stderr, "Error: Energy H_e is NaN for replica %f. AFTER ledn1 SBPmOVE\n", rep->beta);
        //                 exit(1);

        //      }

        RED_n1 = gsl_rng_uniform_int(local_rng, N-rep->boundbp[27]) + rep->boundbp[27];
        perform_sbp_move(RED_n1, 
                        stepsize_a, 
                        stepsize_t, 
                        rep->beta,
                        rep->state,
                        rep->trial,
                        rep->BPSEnergies,
                        rep->BPSEnergiesTrial,
                        rep->BPSEnergies2,
                        rep->BPSEnergiesTrial2,
                        rep->midplanes,
                        rep->mu,
                        &rep->energy_A,
                        &rep->energy_B,
                        &rep->H_e,
                        &sbp_tries,
                        &midbp_tries,
                        &sbp_accept,
                        &midbp_accept, 
                        local_rng, rep->boundbp); 
        // if (isnan(rep->H_e)) {
        //     fprintf(stderr, "Error: Energy H_e is NaN for replica %f. AFTER redn1 SBPmOVE\n", rep->beta);
        //         exit(1);

        //      }
        ///////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////
        ////////////////////////CRANKSHAFT MOVE AND CLUSTER TRANS /////////////
        ///////////////////////////////////////////////////////////////////////
        

        // printf("Starting Crank Moves........................................ \n");
        for (int c=0;c<C_Max;c++) {
            performCrankshaftMoves(stepsize_a, 
                                rep->beta,
                                rep->state,
                                rep->trial,
                                rep->BPSEnergies,
                                rep->BPSEnergiesTrial,
                                rep->BPSEnergies2,
                                rep->BPSEnergiesTrial2,
                                rep->mu, 
                                &rep->energy_A,
                                &rep->energy_B,
                                &rep->H_e,
                                &tries_crank, 
                                &accept_crank, local_rng, rep->boundbp);

            // if (isnan(rep->H_e)) {
            //     fprintf(stderr, "Error: Energy H_e is NaN for replica %f. AFTER Shaft Move\n", rep->beta);
            //         exit(1);

            //      }
                // printf("Starting Crank Trans Moves........................................ \n");

            perform_cluster_trans_moves(stepsize_t, 
                                        rep->beta,
                                        rep->mu,
                                        rep->state,
                                        rep->trial,
                                        rep->BPSEnergies,
                                        rep->BPSEnergiesTrial,
                                        rep->BPSEnergies2,
                                        rep->BPSEnergiesTrial2,
                                        &rep->energy_A,
                                        &rep->energy_B,
                                        &rep->H_e,
                                        &tries_cluster, 
                                        &accept_cluster, local_rng, rep->boundbp);
        }
        

        
        // if (isnan(rep->H_e)) {
        //     fprintf(stderr, "Error: Energy H_e is NaN for replica %f. AFTER Trans Movee\n", rep->beta);
        //         exit(1);

        //      }

    if (signal_was_caught(  )) break;

    
    }

    // printf("Sweep Done........................................ \n");
    // // ONE SWEEP DONE

    // double acceptance_rate_sbp = (double)sbp_accept / (double)sbp_tries;
    // double acceptance_rate_midbp = (double)midbp_accept / (double)midbp_tries;
    // double acceptance_rate_crank = (double)accept_crank / (double)tries_crank;
    // double acceptance_rate_cluster = (double)accept_cluster / (double)tries_cluster;

    // printf("Replica : %f, SBP_RATE: %f, MIDBP_RATE:%f, CRANK_RATE:%f, CLUSTER_RATE:%f\n",
    //  rep->beta, acceptance_rate_sbp, acceptance_rate_midbp, acceptance_rate_crank, acceptance_rate_cluster); 
        

}

















// The swapping function for parallel tempering:
void performReplicaSwaps(Replica **replicas, int num_replicas) {

    int num_swaps = 100;
    int i, j;
    // gsl_rng* local_rg = gsl_rng_alloc(gsl_rng_mt19937);
    double beta_i, beta_j, Ea_i, Eb_i, Ea_j, Eb_j, deltaE_i, deltaE_j, mu_i,mu_j, d, acceptance, delta;

    for (int move = 0; move < num_swaps; move++){
        

         // Randomly choose an adjacent pair (i, i+1)
        i = gsl_rng_uniform_int(rg, num_replicas - 1);  // i in [0, num_replicas-2]
        j = i + 1;  // Adjacent replica index



        // Get beta and energy for replica i and i+1:
        // #pragma omp critical(read_fields)
        // {
        mu_i = replicas[i]->mu;
        mu_j = replicas[j]->mu;
        Ea_i = replicas[i]->energy_A;
        Eb_i = replicas[i]->energy_B;
        Ea_j = replicas[j]->energy_A;
        Eb_j = replicas[j]->energy_B;
        deltaE_i = Ea_i - Eb_i;
        deltaE_j = Ea_j - Eb_j;
        // }
        
        // Compute the swap acceptance probability:
        delta = (mu_i - mu_j) * ( deltaE_i- deltaE_j);
        acceptance = exp(-delta);
        if (acceptance > 1.0) {
            acceptance = 1.0;
        }
        
        // Draw a random number between 0 and 1:
        d = gsl_rng_uniform(rg);
        
        // #pragma omp atomic
        replicas[i]->swap_tries++;
        // #pragma omp atomic
        replicas[i+1]->swap_tries++;

        // If the random number is less than the acceptance probability, swap:
        if (d < acceptance) {
            // Swap the state pointers (i.e., the configurations)

            // #pragma omp critical(replica_swap)
            // {
            DNAState *temp_state = replicas[i]->state;
            replicas[i]->state = replicas[i+1]->state;
            replicas[i+1]->state = temp_state;
            

            DNAState *temp_trial = replicas[i]->trial;
            replicas[i]->trial = replicas[i+1]->trial;
            replicas[i+1]->trial = temp_trial;
        
            
            double *temp_BPSEnergies = replicas[i]->BPSEnergies;
            replicas[i]->BPSEnergies = replicas[i+1]->BPSEnergies;
            replicas[i+1]->BPSEnergies = temp_BPSEnergies;

            double *temp_BPSEnergies2 = replicas[i]->BPSEnergies2;
            replicas[i]->BPSEnergies2 = replicas[i+1]->BPSEnergies2;
            replicas[i+1]->BPSEnergies2 = temp_BPSEnergies2;

            double *temp_BPSEnergiesTrial = replicas[i]->BPSEnergiesTrial;
            replicas[i]->BPSEnergiesTrial = replicas[i+1]->BPSEnergiesTrial;
            replicas[i+1]->BPSEnergiesTrial = temp_BPSEnergiesTrial;

            double *temp_BPSEnergiesTrial2 = replicas[i]->BPSEnergiesTrial2;
            replicas[i]->BPSEnergiesTrial2 = replicas[i+1]->BPSEnergiesTrial2;
            replicas[i+1]->BPSEnergiesTrial2 = temp_BPSEnergiesTrial2;

            RunningStats *temp_energy_stats = &replicas[i]->energy_stats;
            replicas[i]->energy_stats = replicas[i+1]->energy_stats;
            replicas[i+1]->energy_stats = *temp_energy_stats;

            // Swap the energy values:

            double temp_energy_A = replicas[i]->energy_A;
            replicas[i]->energy_A = replicas[i+1]->energy_A;
            replicas[i+1]->energy_A = temp_energy_A;


            double temp_energy_B = replicas[i]->energy_B;
            replicas[i]->energy_B = replicas[i+1]->energy_B;
            replicas[i+1]->energy_B = temp_energy_B;

            // double temp_energy = replicas[i]->H_e;
            // replicas[i]->H_e = replicas[i+1]->H_e;
            // replicas[i+1]->H_e = temp_energy;

            replicas[i]->H_e = (1-replicas[i]->mu)*replicas[i]->energy_A+replicas[i]->mu*replicas[i]->energy_B;
            replicas[i+1]->H_e = (1-replicas[i+1]->mu)*replicas[i+1]->energy_A+replicas[i+1]->mu*replicas[i+1]->energy_B;

            // }




            // Increase the swap counters:
            // #pragma omp atomic
            replicas[i]->accepted_swaps++;
            // #pragma omp atomic
            replicas[i+1]->accepted_swaps++;
            
            // Optional: Print or log that a swap occurred.
            // printf("Swapped replicas %d and %d (acceptance = %f)\n", i, i+1, acceptance);
        }
    }


    // gsl_rng_free(local_rg);
}








// Function to shuffle an integer array using Fisher–Yates shuffle
void shuffle_array(int *array, int n, gsl_rng *r) {
  
    for (int i = n - 1; i > 0; i--) {
        int j = gsl_rng_uniform_int(r, i + 1);
        int temp = array[i];
        array[i] = array[j];
        array[j] = temp;
    }
}






void mc_steps_shuffle(Replica *rep) {
    

    long int sbp_tries = 0; 
    long int midbp_tries = 0;

    long int sbp_accept = 0; 
    long int midbp_accept = 0;

    long int tries_crank = 0;
    long int tries_cluster = 0; 

    long int accept_crank = 0;
    long int accept_cluster = 0; 

    int n1, id_S, id_E, LED_n1, RED_n1;

    double b = rep->beta;
    double stepsize_t = 1.0/(sqrt(b)*10.0);
    double stepsize_a = PI/(sqrt(b)*300.0);
    
    gsl_rng *local_rng = rep->rn;



    long int sweep_size = 20;

    // Define the counts for each move type
    int num_local   = 200;  // local moves
    int num_led     = 1;    // LED move
    int num_red     = 1;    // RED move
    int num_crank   = 10;   // crankshaft moves
    int num_cluster = 10;   // cluster moves
    int move;
    int total_moves = num_local + num_led + num_red + num_crank + num_cluster; 

    // Create an array on the stack for move types
    int moves[total_moves];
    int idx = 0;

    // Fill the array with the appropriate move types
    for (int i = 0; i < num_local; i++) {
        moves[idx++] = MOVE_LOCAL;
    }
    moves[idx++] = MOVE_LED;
    moves[idx++] = MOVE_RED;
    for (int i = 0; i < num_crank; i++) {
        moves[idx++] = MOVE_CRANKSHAFT;
    }
    for (int i = 0; i < num_cluster; i++) {
        moves[idx++] = MOVE_CLUSTER;
    }
 
    
    for (int l = 0; l < sweep_size; l++) {

        // printf("Sweep %d\n", l);
        shuffle_array(moves, total_moves, local_rng);

        // printf("Moves array: ");
        // for (int i = 0; i < total_moves; i++) {
        //     printf("%d ", moves[i]);
        // }
        // printf("\n");

        for (int m = 0; m < total_moves; m++) {

            move = moves[m];

            if (move == MOVE_LOCAL) {
                n1 = gsl_rng_uniform_int(local_rng, N);
                perform_sbp_move(n1, 
                                stepsize_a, 
                                stepsize_t, 
                                rep->beta,
                                rep->state,
                                rep->trial,
                                rep->BPSEnergies,
                                rep->BPSEnergiesTrial,
                                rep->BPSEnergies2,
                                rep->BPSEnergiesTrial2,
                                rep->midplanes,
                                rep->mu,
                                &rep->energy_A,
                                &rep->energy_B,
                                &rep->H_e,
                                &sbp_tries,
                                &midbp_tries,
                                &sbp_accept,
                                &midbp_accept, 
                                local_rng, 
                                rep->boundbp); 

                    
            } else if (move == MOVE_LED) {
               
                LED_n1 = gsl_rng_uniform_int(local_rng, rep->boundbp[0]);
                perform_sbp_move(LED_n1, 
                    stepsize_a, 
                    stepsize_t, 
                    rep->beta,
                    rep->state,
                    rep->trial,
                    rep->BPSEnergies,
                    rep->BPSEnergiesTrial,
                    rep->BPSEnergies2,
                    rep->BPSEnergiesTrial2,
                    rep->midplanes,
                    rep->mu,
                    &rep->energy_A,
                    &rep->energy_B,
                    &rep->H_e,
                    &sbp_tries,
                    &midbp_tries,
                    &sbp_accept,
                    &midbp_accept, local_rng, rep->boundbp);
        
               
        
               
            } else if (move == MOVE_RED) {
             
                RED_n1 = gsl_rng_uniform_int(local_rng, N-rep->boundbp[27]) + rep->boundbp[27];
                perform_sbp_move(RED_n1, 
                                stepsize_a, 
                                stepsize_t, 
                                rep->beta,
                                rep->state,
                                rep->trial,
                                rep->BPSEnergies,
                                rep->BPSEnergiesTrial,
                                rep->BPSEnergies2,
                                rep->BPSEnergiesTrial2,
                                rep->midplanes,
                                rep->mu,
                                &rep->energy_A,
                                &rep->energy_B,
                                &rep->H_e,
                                &sbp_tries,
                                &midbp_tries,
                                &sbp_accept,
                                &midbp_accept, 
                                local_rng, rep->boundbp); 

            } else if (move == MOVE_CRANKSHAFT) {
               performCrankshaftMoves(stepsize_a, 
                                rep->beta,
                                rep->state,
                                rep->trial,
                                rep->BPSEnergies,
                                rep->BPSEnergiesTrial,
                                rep->BPSEnergies2,
                                rep->BPSEnergiesTrial2,
                                rep->mu, 
                                &rep->energy_A,
                                &rep->energy_B,
                                &rep->H_e,
                                &tries_crank, 
                                &accept_crank, local_rng, rep->boundbp);

            } else if (move == MOVE_CLUSTER) {
                
                perform_cluster_trans_moves(stepsize_t, 
                    rep->beta,
                    rep->mu,
                    rep->state,
                    rep->trial,
                    rep->BPSEnergies,
                    rep->BPSEnergiesTrial,
                    rep->BPSEnergies2,
                    rep->BPSEnergiesTrial2,
                    &rep->energy_A,
                    &rep->energy_B,
                    &rep->H_e,
                    &tries_cluster, 
                    &accept_cluster, local_rng, rep->boundbp);
            }
            
            if (signal_was_caught()) break;
        }
        

        
        // if (isnan(rep->H_e)) {
        //     fprintf(stderr, "Error: Energy H_e is NaN for replica %f. AFTER Trans Movee\n", rep->beta);
        //         exit(1);

        //      }

    if (signal_was_caught(  )) break;

    
    }

    // printf("Sweep Done........................................ \n");
    // // ONE SWEEP DONE

    // double acceptance_rate_sbp = (double)sbp_accept / (double)sbp_tries;
    // double acceptance_rate_midbp = (double)midbp_accept / (double)midbp_tries;
    // double acceptance_rate_crank = (double)accept_crank / (double)tries_crank;
    // double acceptance_rate_cluster = (double)accept_cluster / (double)tries_cluster;

    // printf("Replica : %f, SBP_RATE: %f, MIDBP_RATE:%f, CRANK_RATE:%f, CLUSTER_RATE:%f\n",
    //  rep->beta, acceptance_rate_sbp, acceptance_rate_midbp, acceptance_rate_crank, acceptance_rate_cluster); 
        

}









