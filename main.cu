
//Created by Žiga Kos on Jan 18, 2021. Copyright 2021. All rights reserved.
//ziga.kos@fmf.uni-lj.si, zigakos@mit.edu
//Compile with gcc -o prog main.c -lm -lpthread -march=native -std=gnu99 -fopenmp -Wall

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

//Allocation of global matrices
double **U, **F, **FEQ, **FNEW, **P, **Q, **QNEW, **H, **SIGMA;
double *WKONST, *ACTIVITY;
int **E;
char *LMARK;

#include "headers.h"
#include "parameters.c"
#include "utility_functions.c"
#include "lattice_Boltzmann.cu"
#include "finite_difference.cu"


int main(int argc, char** args){
    clock_t start_time, end_time;
    double time_elapsed;
    start_time = clock();

    //Allocation of data arrays
    cudaMallocManaged(&U, 3 * sizeof(double*));  //velocity field; 0th component is density
    for (int m = 0; m < 3; m++) cudaMallocManaged(&U[m], NMAX * sizeof(double));

    cudaMallocManaged(&F, LATTICE_VELOCITY_NUMBER * sizeof(double*));  //Lattice Boltzmann distribution functions
    for (int m = 0; m < LATTICE_VELOCITY_NUMBER; m++) cudaMallocManaged(&F[m], NMAX * sizeof(double));

    cudaMallocManaged(&FNEW, LATTICE_VELOCITY_NUMBER * sizeof(double*));  //matrix field to temporarily store currently calculated F
    for (int m = 0; m < LATTICE_VELOCITY_NUMBER; m++) cudaMallocManaged(&FNEW[m], NMAX * sizeof(double));

    cudaMallocManaged(&FEQ, LATTICE_VELOCITY_NUMBER * sizeof(double*));  //Lattice Boltzmann equilibrium distribution functions
    for (int m = 0; m < LATTICE_VELOCITY_NUMBER; m++) cudaMallocManaged(&FEQ[m], NMAX * sizeof(double));

    cudaMallocManaged(&P, LATTICE_VELOCITY_NUMBER * sizeof(double*));  // Lattice Boltzmann forcing terms
    for (int m = 0; m < LATTICE_VELOCITY_NUMBER; m++) cudaMallocManaged(&P[m], NMAX * sizeof(double));

    cudaMallocManaged(&Q, 2 * sizeof(double*));  //Q-tensor
    for (int m = 0; m < 2; m++) cudaMallocManaged(&Q[m], NMAX * sizeof(double));

    cudaMallocManaged(&QNEW, 2 * sizeof(double*));  //Q-tensor to temporarily store values
    for (int m = 0; m < 2; m++) cudaMallocManaged(&QNEW[m], NMAX * sizeof(double));

    cudaMallocManaged(&H, 2 * sizeof(double*));  //Molecular field (from Q-tensor relaxation)
    for (int m = 0; m < 2; m++) cudaMallocManaged(&H[m], NMAX * sizeof(double));

    cudaMallocManaged(&SIGMA, 4 * sizeof(double*));  //Stress tensor
    for (int m = 0; m < 4; m++) cudaMallocManaged(&SIGMA[m], NMAX * sizeof(double));

    cudaMallocManaged(&E, LATTICE_VELOCITY_NUMBER * sizeof(int*));  // matrix of all LB characteristic velocity lattice vectors
    for (int m = 0; m < LATTICE_VELOCITY_NUMBER; m++) cudaMallocManaged(&E[m], 2 * sizeof(int));

    cudaMallocManaged(&WKONST, LATTICE_VELOCITY_NUMBER * sizeof(double));  //Lattice Boltzmann weights
    cudaMallocManaged(&LMARK, NMAX * sizeof(char));        //Logical markers to determine a mesh point function (bulk or boundary condition)
    cudaMallocManaged(&ACTIVITY, NMAX * sizeof(double));  // activity vs. no activity

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

    for (int l = 0; l < NMAX; l++) {
        //Logical markers
        if ((l % I) == 0 || (l % I) == I - 1 || ((l % (I * J)) / I) == 0 || ((l % (I * J)) / I) == J - 1) LMARK[l] = LMARKBC;
        else LMARK[l] = LMARKBULK;

        if (isPointInActivityPattern(l)) ACTIVITY[l] = ALPHA;
        else ACTIVITY[l] = 0.0;
        
        //Velocity Field
        double angle = 0.01 * (double)rand() / (double)((unsigned)RAND_MAX + 1);  // randomly initialize velocity field
        U[0][l] = DENSITYINIT;
        U[1][l] = INLET_VELOCITY; //0.0; //0.001 * cos(angle);
        U[2][l] = 0.0; //0.001 * sin(angle);
        
        //Q tensor
        double dx_defect = (double)(l % I) - defect_x;
        double dy_defect = (double)((l % (I * J)) / I) - defect_y;
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
    int num_blocks = (NMAX + BLOCK_SIZE - 1) / BLOCK_SIZE;  // for CUDA

    computeFeq<<<num_blocks, BLOCK_SIZE>>>(U, E, FEQ, WKONST);
    cudaDeviceSynchronize();

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
        compute_FD_step(U, Q, QNEW, H, LMARK);
        compute_LB_step(U, E, FNEW, F, FEQ, P, Q, H, SIGMA, WKONST, ACTIVITY, LMARK);
    }
    
    //Write the final fields
    write_velocity(TIME_STEPS);
    write_orientation(TIME_STEPS);
    
    
    //Free up space
    for (int m = 0; m < 3; m++) cudaFree(U[m]);
    cudaFree(U);

    for (int m = 0; m < LATTICE_VELOCITY_NUMBER; m++) cudaFree(F[m]);
    cudaFree(F);

    for (int m = 0; m < LATTICE_VELOCITY_NUMBER; m++) cudaFree(FNEW[m]);
    cudaFree(FNEW);

    for (int m = 0; m < LATTICE_VELOCITY_NUMBER; m++) cudaFree(FEQ[m]);
    cudaFree(FEQ);

    for (int m = 0; m < LATTICE_VELOCITY_NUMBER; m++) cudaFree(P[m]);
    cudaFree(P);

    for (int m = 0; m < 2; m++) cudaFree(Q[m]);
    cudaFree(Q);

    for (int m = 0; m < 2; m++) cudaFree(QNEW[m]);
    cudaFree(QNEW);

    for (int m = 0; m < 2; m++) cudaFree(H[m]);
    cudaFree(H);

    for (int m = 0; m < 4; m++) cudaFree(SIGMA[m]);
    cudaFree(SIGMA);

    cudaFree(WKONST);

    for (int m = 0; m < LATTICE_VELOCITY_NUMBER; m++) cudaFree(E[m]);
    cudaFree(E);

    cudaFree(LMARK);

    cudaFree(ACTIVITY);

    end_time = clock();
    time_elapsed = ((double) (end_time - start_time)) / CLOCKS_PER_SEC;
    printf("Simulation took %f seconds to complete\n", time_elapsed);

    return 0;
}
