
//Created by Žiga Kos on Jan 18, 2021. Copyright 2021. All rights reserved.
//ziga.kos@fmf.uni-lj.si, zigakos@mit.edu
//Compile with gcc -o prog main.c -lm -lpthread -march=native -std=gnu99 -fopenmp -Wall

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>

//Allocation of global matrices
double **U, **F, **FEQ, **FNEW, **P;
double **Q, **QNEW;
double *WKONST;
int **E;
char *LMARK;
double *ACTIVITY;
double **H;
double **SIGMA;

#include "headers.h"
#include "parameters.c"
#include "utility_functions.c"
#include "lattice_Boltzmann.c"
#include "finite_difference.c"



int main(int argc, char** args){
    
    //Allocation of data arrays
    U = (double**)malloc(3 * sizeof(double*));  //velocity field; 4th component is density
    for (int m = 0; m < 3; m++) U[m] = (double*)malloc(NMAX * sizeof(double));

    F = (double**)malloc(LATTICE_VELOCITY_NUMBER * sizeof(double*));  //Lattice Boltzmann distribution functions
    for (int m = 0; m < LATTICE_VELOCITY_NUMBER; m++) F[m] = (double*)malloc(NMAX * sizeof(double));

    FNEW = (double**)malloc(LATTICE_VELOCITY_NUMBER * sizeof(double*));  //matrix field to temporarily store currently calculated F
    for (int m = 0; m < LATTICE_VELOCITY_NUMBER; m++) FNEW[m] = (double*)malloc(NMAX * sizeof(double));

    FEQ = (double**)malloc(LATTICE_VELOCITY_NUMBER * sizeof(double*));  //Lattice Boltzmann equilibrium distribution functions
    for (int m = 0; m < LATTICE_VELOCITY_NUMBER; m++) FEQ[m] = (double*)malloc(NMAX * sizeof(double));

    P = (double**)malloc(LATTICE_VELOCITY_NUMBER * sizeof(double*));  // Lattice Boltzmann forcing terms
    for (int m = 0; m < LATTICE_VELOCITY_NUMBER; m++) P[m] = (double*)malloc(NMAX * sizeof(double));

    Q = (double**)malloc(2 * sizeof(double*));  //Q-tensor
    for (int m = 0; m < 2; m++) Q[m] = (double*)malloc(NMAX * sizeof(double));

    QNEW = (double**)malloc(2 * sizeof(double*));  //Q-tensor to temporarily store values
    for (int m = 0; m < 2; m++) QNEW[m] = (double*)malloc(NMAX * sizeof(double));

    H = (double**)malloc(2 * sizeof(double*));  //Molecular field (from Q-tensor relaxation)
    for (int m = 0; m < 2; m++) H[m] = (double*)malloc(NMAX * sizeof(double));

    SIGMA = (double**)malloc(4 * sizeof(double*));  //Stress tensor
    for (int m = 0; m < 4; m++) SIGMA[m] = (double*)malloc(NMAX * sizeof(double));

    E = (int**)malloc(LATTICE_VELOCITY_NUMBER * sizeof(int*));  // matrix of all LB characteristic velocity lattice vectors
    for (int m = 0; m < LATTICE_VELOCITY_NUMBER; m++) E[m] = (int*)malloc(2 * sizeof(int));

    WKONST = (double*)malloc(LATTICE_VELOCITY_NUMBER * sizeof(double));  //Lattice Boltzmann weights
    LMARK = malloc(NMAX * sizeof(char));        //Logical markers to determine a mesh point function (bulk or boundary condition)
    ACTIVITY = malloc(NMAX * sizeof(double));  // activity vs. no activity

    //Initialization
    unsigned seed = 12345;
    srand(seed);
    initialiseE(); //Initializes the lattice vectors

    // defect location
    double defect_x = round(0.05 * I + 0.5);
    //double defect_x = round(0.24 * I - 0.5);
    double defect_y = round(0.5 * J - 0.5);
    double phi0 = 0.5 * M_PI;   // 0.5 * M_PI;
    double topo_charge = 0.5;

    int obs_top_row = round(0.25 * J - 0.5);
    int obs_bot_row = round(0.25 * J - 0.5);
    int obs_left_col = round(0.20 * I - 0.5);
    int obs_right_col = round(0.40 * I - 0.5);

    for (int l = 0; l < NMAX; l++) {
        //Logical markers
        if (j_vr(l) < obs_top_row && j_vr(l) > 0 && i_vr(l) == obs_left_col) {
            LMARK[l] = LMARKOBSLEFT;
        } else if (j_vr(l) < obs_top_row && j_vr(l) > 0 && i_vr(l) == obs_right_col) {
            LMARK[l] = LMARKOBSRIGHT;
        } else if (j_vr(l) < obs_top_row && i_vr(l) > obs_left_col && i_vr(l) < obs_right_col) {
            LMARK[l] = LMARKOBSBULK;
        } else if (j_vr(l) == obs_top_row && i_vr(l) > obs_left_col && i_vr(l) < obs_right_col) {
            LMARK[l] = LMARKOBSTOP;
        } else if (j_vr(l) == obs_top_row && i_vr(l) == obs_left_col) {
            LMARK[l] = LMARKOBS_TOP_LEFT_CORNER;
        } else if (j_vr(l) == obs_top_row && i_vr(l) == obs_right_col) {
            LMARK[l] = LMARKOBS_TOP_RIGHT_CORNER;
        } else if (j_vr(l) == 0 && i_vr(l) == obs_left_col) {
            LMARK[l] = LMARKOBS_BOT_LEFT_CORNER;
        } else if (j_vr(l) == 0 && i_vr(l) == obs_right_col) {
            LMARK[l] = LMARKOBS_BOT_RIGHT_CORNER;
        } else if (i_vr(l) == 0 || i_vr(l) == I - 1 || j_vr(l) == 0 || j_vr(l) == J - 1) {
            LMARK[l] = LMARKBC;
        } else {
            LMARK[l] = LMARKBULK;
        }

        if (isPointInActivityPattern(l)) ACTIVITY[l] = ALPHA;
        else ACTIVITY[l] = 0.0;
        
        //Velocity Field
        double angle = 0.01 * (double)rand() / (double)((unsigned)RAND_MAX + 1);  // randomly initialize velocity field
        U[0][l] = DENSITYINIT;
        U[1][l] = INLET_VELOCITY; //0.001 * cos(angle);
        U[2][l] = 0.0; //0.001 * sin(angle);
        
        //Q tensor
        double dx_defect = (double)i_vr(l) - defect_x;
        double dy_defect = (double)j_vr(l) - defect_y;
        //angle = phi0 + topo_charge * atan2(dy_defect, dx_defect);
        angle = phi0;

        if (dx_defect <= 0) {
            angle += (1.0 - abs(dy_defect) / defect_y) * topo_charge * atan2(dy_defect, dx_defect);  // TODO: generalize this to work with the defect not being at the y midpoint
        }
        //angle = phi0 + cos(0.5 * M_PI * abs(dy_defect) / defect_y) * topo_charge * atan2(dy_defect, dx_defect);
        //angle = M_PI * (double)rand() / (double)((unsigned)RAND_MAX + 1);  // randomly initialize Q tensor

        // NOTE:
        // Q = s * ([[ cos^2(theta) - 1/2,      cos(theta) * sin(theta) ],
        //           [ cos(theta) * sin(theta), sin^2(theta) - 1/2      ]])
        //   = (s/2) * ([[ cos(2 * theta),  sin(2 * theta)]],
        //               [ sin(2 * theta), -cos(2 * theta)]])      // using trig identities
        //
        // We only store two values for Q (Q11 and Q12) since Q22 = -Q11 and Q12 = Q21

        double degree_of_order = 1.;  // degree of order (s) must be between -1/2 and 1
        Q[0][l] = degree_of_order / 2.0 * cos(2 * angle);  //Qxx component
        Q[1][l] = degree_of_order / 2.0 * sin(2 * angle);  //Qxy component
    }
    
    //Compute distribution functions from velocity initialization
    computeFeq();
    #pragma omp parallel for num_threads(STPROC) schedule(dynamic)
    for (int l = 0; l < NMAX; l++) {
        for(int m = 0; m < LATTICE_VELOCITY_NUMBER; m++) {
            F[m][l] = FEQ[m][l];
            FNEW[m][l] = F[m][l];
        }
    }
    
    //Main part - time evolution
    for (int t = 0; t < TIME_STEPS; t++) {
        if (t % TIME_WRITE == 0) {
            write_velocity(t);
            write_orientation(t);
        }
        compute_FD_step();
        compute_LB_step();
    }
    
    //Write the final fields
    write_velocity(TIME_STEPS);
    write_orientation(TIME_STEPS);
    
    
    //Free up space
    for (int m = 0; m < 3; m++) free(U[m]);
    free(U);

    for (int m = 0; m < LATTICE_VELOCITY_NUMBER; m++) free(F[m]);
    free(F);

    for (int m = 0; m < LATTICE_VELOCITY_NUMBER; m++) free(FNEW[m]);
    free(FNEW);

    for (int m = 0; m < LATTICE_VELOCITY_NUMBER; m++) free(FEQ[m]);
    free(FEQ);

    for (int m = 0; m < LATTICE_VELOCITY_NUMBER; m++) free(P[m]);
    free(P);

    for (int m = 0; m < 2; m++) free(Q[m]);
    free(Q);

    for (int m = 0; m < 2; m++) free(QNEW[m]);
    free(QNEW);

    for (int m = 0; m < 2; m++) free(H[m]);
    free(H);

    for (int m = 0; m < 4; m++) free(SIGMA[m]);
    free(SIGMA);

    free(WKONST);

    for (int m = 0; m < LATTICE_VELOCITY_NUMBER; m++) free(E[m]);
    free(E);

    free(LMARK);

    free(ACTIVITY);
    
    return 0;
}
