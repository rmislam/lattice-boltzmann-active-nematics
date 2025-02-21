//Makes a main step in the LB part
void compute_LB_step(double **U, int **E, double **FNEW, double **F, double **FEQ, double **P, double **Q, double **H, double **SIGMA, double *WKONST, double *ACTIVITY, char *LMARK) {
    int num_blocks = (NMAX + BLOCK_SIZE - 1) / BLOCK_SIZE;

    computeFeq<<<num_blocks, BLOCK_SIZE>>>(U, E, FEQ, WKONST);   //Compute the equilibrium distribution

    computeSigma<<<num_blocks, BLOCK_SIZE>>>(SIGMA, Q, H, ACTIVITY, LMARK);  //Compute stress tensor
    cudaDeviceSynchronize();

    computeP<<<num_blocks, BLOCK_SIZE>>>(U, E, P, SIGMA, WKONST, LMARK);  //Compute the forcing terms
    cudaDeviceSynchronize();

    computeFNEW<<<num_blocks, BLOCK_SIZE>>>(FNEW, F, FEQ, E, P);  // Compute the new distribution values (i.e., apply the Boltzmann equation)
    cudaDeviceSynchronize();

    enforceBCAndCalcFNEW2F2U<<<num_blocks, BLOCK_SIZE>>>(FNEW, F, U, E, SIGMA, LMARK);
    cudaDeviceSynchronize();
}

//Computes equilibrium distribution function
__global__
void computeFeq(double **U, int **E, double **FEQ, double *WKONST) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;

    for (int l = index; l < NMAX; l += stride) {
        double u2 = U[1][l] * U[1][l] + U[2][l] * U[2][l];	//Velocity squared

        for (int m = 0; m < LATTICE_VELOCITY_NUMBER; m++) {
            double ue = U[1][l] * E[m][0] + U[2][l] * E[m][1];		//Product of velocity and characteristic vectors
            FEQ[m][l] = WKONST[m] * U[0][l] * (1.0 + 3.0 * ue + 9.0 / 2.0 * ue * ue - 1.5 * u2);  // NOTE: no need to modify
        }
    }
}

//Compute stress tensor
__global__
void computeSigma(double **SIGMA, double **Q, double **H, double *ACTIVITY, char *LMARK) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;

    for (int l = index; l < NMAX; l += stride) {
        if (!(LMARK[l] == LMARK_BULK)) continue;
        // NOTE: positive activity is extensile, negative activity is contractile
        // NOTE: Qxx = Q[0], Qyy = -Q[0], Qxy = Qyx = Q[1], and Hxx = H[0], Hyy = -H[0], Hxy = Hyx = H[1]
        double dxQ0 = (Q[0][l + 1] - Q[0][l - 1]) / 2.0;
        double dxQ1 = (Q[1][l + 1] - Q[1][l - 1]) / 2.0;
        double dyQ0 = (Q[0][l + I] - Q[0][l - I]) / 2.0;
        double dyQ1 = (Q[1][l + I] - Q[1][l - I]) / 2.0;
        double gradQsq = dxQ0 * dxQ0 + dxQ1 * dxQ1 + dyQ0 * dyQ0 + dyQ1 * dyQ1;

        //xx component of the stress tensor
        SIGMA[0][l] = L * gradQsq + 4 * XI * (Q[0][l] * H[0][l] + Q[1][l] * H[1][l]) * (Q[0][l] + 0.5) - XI * (2 * (Q[0][l] * H[0][l] + Q[1][l] * H[1][l]) + H[0][l]) - 2 * L * (dxQ0 * dxQ0 + dxQ1 * dxQ1);
        //xy component of the stress tensor
        SIGMA[1][l] = 4 * XI * (Q[0][l] * H[0][l] + Q[1][l] * H[1][l]) * Q[1][l] - XI * H[1][l] - 2 * L * (dxQ0 * dyQ0 + dxQ1 * dyQ1) + 2 * (Q[0][l] * H[1][l] - Q[1][l] * H[0][l]);
        //yx component of the stress tensor
        SIGMA[2][l] = 4 * XI * (Q[0][l] * H[0][l] + Q[1][l] * H[1][l]) * Q[1][l] - XI * H[1][l] - 2 * L * (dxQ0 * dyQ0 + dxQ1 * dyQ1) + 2 * (Q[1][l] * H[0][l] - Q[0][l] * H[1][l]);
        //yy component of the stress tensor
        SIGMA[3][l] = L * gradQsq + 4 * XI * (Q[0][l] * H[0][l] + Q[1][l] * H[1][l]) * (-Q[0][l] + 0.5) - XI * (2 * (Q[0][l] * H[0][l] + Q[1][l] * H[1][l]) - H[0][l]) - 2 * L * (dyQ0 * dyQ0 + dyQ1 * dyQ1);

        // add activity
        SIGMA[0][l] -= ACTIVITY[l] * Q[0][l];
        SIGMA[1][l] -= ACTIVITY[l] * Q[1][l];
        SIGMA[2][l] -= ACTIVITY[l] * Q[1][l];
        SIGMA[3][l] += ACTIVITY[l] * Q[0][l];
    }
}

//Compute the forcing terms
__global__
void computeP(double **U, int **E, double **P, double **SIGMA, double *WKONST, char *LMARK) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;

    for (int l = index; l < NMAX; l += stride) {
        if (!(LMARK[l] == LMARK_BULK)) continue;
        //Compute the derivatives of the stress tensor and the force
        double forceX = (SIGMA[0][l + 1] - SIGMA[0][l - 1]) / 2.0 + (SIGMA[1][l + I] - SIGMA[1][l - I]) / 2.0 - MU * U[1][l];
        double forceY = (SIGMA[2][l + 1] - SIGMA[2][l - 1]) / 2.0 + (SIGMA[3][l + I] - SIGMA[3][l - I]) / 2.0 - MU * U[2][l];

        double uF = U[1][l] * forceX + U[2][l] * forceY;  //Product of force and velocity
        for (int m = 0; m < LATTICE_VELOCITY_NUMBER; m++) {
            double ue = U[1][l] * E[m][0] + U[2][l] * E[m][1];	//Product of velocity and characteristic vectors
            double eF = E[m][0] * forceX + E[m][1] * forceY;  //Product of force and characteristic vectors
            P[m][l] = (1.0 - DT / 2.0 / TAUF) * WKONST[m] * (3.0 * eF - 3.0 * uF + 9.0 * ue * eF);  // "source term" S_i from Kruger 6.2
        }
    }
}

__global__
void computeFNEW(double **FNEW, double **F, double **FEQ, int **E, double **P) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;

    for (int l = index; l < NMAX; l += stride) {
        for (int m = 0; m < LATTICE_VELOCITY_NUMBER; m++) {
            // Compute streaming location
            int lnew = (l % I) + E[m][0] + (((l % (I * J)) / I) + E[m][1]) * I;  // i + j * I
            if (lnew >= NMAX || lnew < 0) continue;
            FNEW[m][lnew] = F[m][l] + DT * ((FEQ[m][l] - F[m][l]) / TAUF + P[m][l]);
        }
    }
}

__global__
void enforceBCAndCalcFNEW2F2U(double **FNEW, double **F, double **U, int **E, double **SIGMA, char *LMARK) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;

    for (int l = index; l < NMAX; l += stride) {
        // enforce boundary conditions
        if (LMARK[l] == LMARK_TOP_WALL) {  // walls -- no-slip
            FNEW[4][l] = F[2][l];
            FNEW[7][l] = F[5][l];
            FNEW[8][l] = F[6][l];
        } else if (LMARK[l] == LMARK_BOT_WALL) {
            FNEW[2][l] = F[4][l];
            FNEW[5][l] = F[7][l];
            FNEW[6][l] = F[8][l];
        } else if (LMARK[l] == LMARK_RIGHT_WALL) {
            FNEW[3][l] = F[1][l];
            FNEW[6][l] = F[8][l];
            FNEW[7][l] = F[5][l];
        } else if (LMARK[l] == LMARK_LEFT_WALL) {
            FNEW[1][l] = F[3][l];
            FNEW[8][l] = F[6][l];
            FNEW[5][l] = F[7][l];
        } else if (LMARK[l] == LMARK_CORNER_TOP_LEFT) {  // corners -- no-slip
            FNEW[1][l] = F[3][l];
            FNEW[4][l] = F[2][l];
            FNEW[8][l] = F[6][l];
        } else if (LMARK[l] == LMARK_CORNER_BOT_LEFT) {
            FNEW[1][l] = F[3][l];
            FNEW[2][l] = F[4][l];
            FNEW[5][l] = F[7][l];
        } else if (LMARK[l] == LMARK_CORNER_TOP_RIGHT) {
            FNEW[3][l] = F[1][l];
            FNEW[4][l] = F[2][l];
            FNEW[7][l] = F[5][l];
        } else if (LMARK[l] == LMARK_CORNER_BOT_RIGHT) {
            FNEW[3][l] = F[1][l];
            FNEW[2][l] = F[4][l];
            FNEW[6][l] = F[8][l];
        } else if (LMARK[l] == LMARK_TOP_OUTLET) {  // outlets -- periodic
            FNEW[4][l] = F[4][l - I * (J - 2)];
            FNEW[7][l] = F[7][l - I * (J - 2)];
            FNEW[8][l] = F[8][l - I * (J - 2)];
        } else if (LMARK[l] == LMARK_BOT_OUTLET) {
            FNEW[2][l] = F[2][l + I * (J - 2)];
            FNEW[5][l] = F[5][l + I * (J - 2)];
            FNEW[6][l] = F[6][l + I * (J - 2)];
        } else if (LMARK[l] == LMARK_LEFT_OUTLET) {
            FNEW[1][l] = F[1][l + I - 2];
            FNEW[5][l] = F[5][l + I - 2];
            FNEW[8][l] = F[8][l + I - 2];
        } else if (LMARK[l] == LMARK_RIGHT_OUTLET) {
            FNEW[3][l] = F[3][l - (I - 2)];
            FNEW[6][l] = F[6][l - (I - 2)];
            FNEW[7][l] = F[7][l - (I - 2)];
        }

        // Copies FNEW to F
        for (int m = 0; m < LATTICE_VELOCITY_NUMBER; m++) {
            F[m][l] = FNEW[m][l];
        }

        if (LMARK[l] == LMARK_BULK) {
            double fex = 0.0, fey = 0.0, density = 0.0;

            for (int m = 0; m < LATTICE_VELOCITY_NUMBER; m++) {
                density += F[m][l];
                fex += F[m][l] * E[m][0];
                fey += F[m][l] * E[m][1];
            }

            double forceX = (SIGMA[0][l + 1] - SIGMA[0][l - 1]) / 2.0 + (SIGMA[1][l + I] - SIGMA[1][l - I]) / 2.0 - MU * U[1][l];
            double forceY = (SIGMA[2][l + 1] - SIGMA[2][l - 1]) / 2.0 + (SIGMA[3][l + I] - SIGMA[3][l - I]) / 2.0 - MU * U[2][l];

            U[1][l] = fex / density + forceX * DT / 2.0 / density;
            U[2][l] = fey / density + forceY * DT / 2.0 / density;
            U[0][l] = density;
        } else if (LMARK[l] == LMARK_TOP_OUTLET) {  // outlets -- periodic
            U[1][l] = U[1][l - I * (J - 2)];
            U[2][l] = U[2][l - I * (J - 2)];
        } else if (LMARK[l] == LMARK_BOT_OUTLET) {
            U[1][l] = U[1][l + I * (J - 2)];
            U[2][l] = U[2][l + I * (J - 2)];
        } else if (LMARK[l] == LMARK_LEFT_OUTLET) {
            U[1][l] = U[1][l + I - 2];
            U[2][l] = U[2][l + I - 2];
        } else if (LMARK[l] == LMARK_RIGHT_OUTLET) {
            U[1][l] = U[1][l - (I - 2)];
            U[2][l] = U[2][l - (I - 2)];
        } else if (LMARK[l] == LMARK_TOP_WALL || LMARK[l] == LMARK_BOT_WALL || LMARK[l] == LMARK_LEFT_WALL || LMARK[l] == LMARK_RIGHT_WALL) {  // no-slip condition along channel walls
            U[1][l] = 0;
            U[2][l] = 0;
        } else if (LMARK[l] == LMARK_CORNER_BOT_LEFT || LMARK[l] == LMARK_CORNER_BOT_RIGHT || LMARK[l] == LMARK_CORNER_TOP_LEFT || LMARK[l] == LMARK_CORNER_TOP_RIGHT) {  // no-slip on corners
            U[1][l] = 0;
            U[2][l] = 0;
        } else if (LMARK[l] == LMARK_OBS_BULK) {  // nothing moving inside obstacle
            U[1][l] = 0;
            U[2][l] = 0;
        }
    }
}