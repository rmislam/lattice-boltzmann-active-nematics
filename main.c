// Modified and extended by Russ Islam
// Original tutorial code created by Žiga Kos and Miha Ravnik: https://zenodo.org/records/4737814

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

    // channel and obstacle structure
    int bot_row = round(0.4 * J - 0.5);
    int top_row = round(0.6 * J - 0.5);
    int left_col = round(0.4 * I - 0.5);
    int right_col = round(0.6 * I - 0.5);

    for (int l = 0; l < NMAX; l++) {
        int i = i_vr(l);
        int j = j_vr(l);

        //Logical markers
        if (i == 0 && j < top_row && j > bot_row) { // outlets
            LMARK[l] = LMARK_LEFT_OUTLET;
        } else if (i == I - 1 && j < top_row && j > bot_row) {
            LMARK[l] = LMARK_RIGHT_OUTLET;
        } else if (i > left_col && i < right_col && j == J - 1) {
            LMARK[l] = LMARK_TOP_OUTLET;
        } else if (i > left_col && i < right_col && j == 0) {
            LMARK[l] = LMARK_BOT_OUTLET;
        } else if (j == bot_row && (i < left_col || i > right_col)) { // walls
            LMARK[l] = LMARK_BOT_WALL;
        } else if (j == top_row && (i < left_col || i > right_col)) {
            LMARK[l] = LMARK_TOP_WALL;
        } else if (i == left_col && (j < bot_row || j > top_row)) {
            LMARK[l] = LMARK_LEFT_WALL;
        } else if (i == right_col && (j < bot_row || j > top_row)) {
            LMARK[l] = LMARK_RIGHT_WALL;
        } else if (i == left_col && j == bot_row) { // corners
            LMARK[l] = LMARK_CORNER_BOT_LEFT;
        } else if (i == right_col && j == bot_row) {
            LMARK[l] = LMARK_CORNER_BOT_RIGHT;
        } else if (i == left_col && j == top_row) {
            LMARK[l] = LMARK_CORNER_TOP_LEFT;
        } else if (i == right_col && j == top_row) {
            LMARK[l] = LMARK_CORNER_TOP_RIGHT;
        } else if ((i < left_col && j < bot_row) || (i > right_col && j < bot_row) || (i < left_col && j > top_row) || (i > right_col && j > top_row)) { // bulk
            LMARK[l] = LMARK_OBS_BULK;
        } else {
            LMARK[l] = LMARK_BULK;
        }

        //if (isPointInActivityPattern(l)) ACTIVITY[l] = ALPHA;
        //else ACTIVITY[l] = 0.0;
        ACTIVITY[l] = 0.0;  // initially set activity off everywhere. turn it on below
        
        //Velocity Field
        //double angle = 0.01 * (double)rand() / (double)((unsigned)RAND_MAX + 1);  // randomly initialize velocity field
        U[0][l] = DENSITYINIT;
        U[1][l] = INLET_VELOCITY; //0.0; //0.001 * cos(angle);
        U[2][l] = 0.0; //0.001 * sin(angle);

        if (i < left_col || i > right_col) {
            double angle = 0.5 * M_PI;
            Q[0][l] = 1.0 / 2.0 * cos(2 * angle);  //Qxx component
            Q[1][l] = 1.0 / 2.0 * sin(2 * angle);  //Qxy component
        } else if (j < bot_row || j > top_row) {
            double angle = 0;
            Q[0][l] = 1.0 / 2.0 * cos(2 * angle);  //Qxx component
            Q[1][l] = 1.0 / 2.0 * sin(2 * angle);  //Qxy component
        } else {
            // Initialize director to point between quadrants II and IV
            double angle = 0.75 * M_PI;
            Q[0][l] = 1.0 / 2.0 * cos(2 * angle);  //Qxx component
            Q[1][l] = 1.0 / 2.0 * sin(2 * angle);  //Qxy component
        }

        //double angle = 0.75 * M_PI;
        //Q[0][l] = 1.0 / 2.0 * cos(2 * angle);  //Qxx component
        //Q[1][l] = 1.0 / 2.0 * sin(2 * angle);  //Qxy component
    }

    // defect location
    double defect_x = round(0.1 * I + 0.5); //round(0.5 * I - 0.5);
    //double defect_x = round(0.24 * I - 0.5);
    double defect_y = round(0.5 * J - 0.5); //round(0.1 * J + 0.5); 
    double phi0 = 0.5 * M_PI;   // 0.5 * M_PI;
    double topo_charge = 0.5;
    double degree_of_order = 1.;

    //createDefect(defect_x, defect_y, phi0, topo_charge, degree_of_order);
    
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
        if (t == WARM_UP_STEPS) {
            for (int l = 0; l < NMAX; l++) {
                if (isPointInActivityPattern(l)) ACTIVITY[l] = ALPHA;
                else ACTIVITY[l] = 0.0;
            }

            createDefect(defect_x, defect_y, phi0, topo_charge, degree_of_order);
        }

        if (t >= TIME_START_WRITE && t % TIME_WRITE == 0) {
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
