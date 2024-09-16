 //Computes equilibrium distribution function
void computeFeq() {
    #pragma omp parallel for num_threads(STPROC) schedule(dynamic)
    for (int l = 0; l < NMAX; l++) {
        if (!(LMARK[l] == LMARKBULK)) continue;
    
        double u2 = U[1][l] * U[1][l] + U[2][l] * U[2][l];	//Velocity squared
    
        for (int m = 0; m < LATTICE_VELOCITY_NUMBER; m++) {
            double ue = U[1][l] * E[m][0] + U[2][l] * E[m][1];		//Product of velocity and characteristic vectors
            FEQ[m][l] = WKONST[m] * U[0][l] * (1.0 + 3.0 * ue + 9.0 / 2.0 * ue * ue - 1.5 * u2);
        }
    }
}


//Makes a main step in the LB part
void compute_LB_step() {
    computeFeq();   //Compute the equilibrium distribution
    computeP(); //Compute the forcing terms
    computeOBC_LB();    //Makes the open boundary conditions for LB functions
    
    #pragma omp parallel for num_threads(STPROC) schedule(dynamic)
    for (int l = 0; l < NMAX; l++) {
        for (int m = 0; m < LATTICE_VELOCITY_NUMBER; m++) {
            int lnew = calcLBlnew(l, m);   //Streaming location
            if (lnew >= NMAX || lnew < 0) continue;
            FNEW[m][lnew] = F[m][l] - DT / TAUF * (F[m][l] - FEQ[m][l]) + P[m][l];    //Streaming and collision
        }
    }
    
    calcFNEW2F();   //Copies FNEW to F
    
    //Compute the velocity field
    #pragma omp parallel for num_threads(STPROC) schedule(dynamic)
    for (int l = 0; l < NMAX; l++) calcF2U(l);
    
    //Adapt the velocity field for OBC
    #pragma omp parallel for num_threads(STPROC) schedule(dynamic)

    for(int l = 0; l < NMAX; l++) {
        for (int m = 0; m < 3; m++) {
            U[m][l] = U[m][calcLobc(l)];
        }
    }
}


//Compute the forcing terms
void computeP() {
    //Compute stress tensor
    compute_sigma();
    
    #pragma omp parallel for num_threads(STPROC) schedule(dynamic)
    for (int l = 0; l < NMAX; l++) {
        if (!(LMARK[l] == LMARKBULK)) continue;
        //Compute the derivatives of the stress tensor and the force
        double forceX = (SIGMA[0][l + 1] - SIGMA[0][l - 1]) / 2.0 + (SIGMA[1][l + I] - SIGMA[1][l - I]) / 2.0 - MU * U[1][l];
        double forceY = (SIGMA[2][l + 1] - SIGMA[2][l - 1]) / 2.0 + (SIGMA[3][l + I] - SIGMA[3][l - I]) / 2.0 - MU * U[2][l];
        
//        if (i_vr(l) == I/2 && j_vr(l) == 20) printf("%lf %lf    ", forceX, forceY);
//        if (i_vr(l) == I/2 && j_vr(l) == 20) printf("%lf %lf %f %f      ", SIGMA[0][l], SIGMA[1][l], SIGMA[2][l], SIGMA[3][l]);
        
        double uF = U[1][l] * forceX + U[2][l] * forceY;  //Product of force and velocity
        for (int m = 0; m < LATTICE_VELOCITY_NUMBER; m++) {
            double ue = U[1][l] * E[m][0] + U[2][l] * E[m][1];	//Product of velocity and characteristic vectors
            double eF = E[m][0] * forceX + E[m][1] * forceY;  //Product of force and characteristic vectors
            P[m][l] = (1.0 - DT / 2.0 / TAUF) * WKONST[m] * ( 3.0 * eF - 3.0 * uF + 9.0 * ue * eF );
//            P[m][l] = 0;
        }
    }
}


//Enforce open boundaries
void computeOBC_LB() {
    #pragma omp parallel for num_threads(STPROC) schedule(dynamic)
    for (int l = 0; l < NMAX; l++) {
        if (LMARK[l] & LMARKBC) {
            int lobc = calcLobc(l);
            for (int m = 0; m < LATTICE_VELOCITY_NUMBER; m++) {
                P[m][l] = P[m][lobc];
                FEQ[m][l] = FEQ[m][lobc];
                F[m][l] = F[m][lobc];
            }
        }
    }
}

//Coumpute the index of the open boundary
int calcLobc(int l) {
    int xp2 = 0, yp2 = 0, zp2 = 0;

    // open boundaries on left and right, periodic boundaries on top and bottom
    if (i_vr(l) == 0)      xp2 = 1;
    if (i_vr(l) == I - 1)  xp2 = -1;
    if (j_vr(l) == 0)      yp2 = I * (J - 2);
    if (j_vr(l) == J - 1)  yp2 = -I * (J - 2);

    return (l + xp2 + yp2 + zp2);
}

//Computes streaming location
int calcLBlnew(int l, int m) {
    int i, j;
    i = i_vr(l) + E[m][0];
    j = j_vr(l) + E[m][1];
    return i + j * I;
}

//Copies FNEW to F
void calcFNEW2F() {
    #pragma omp parallel for num_threads(STPROC) schedule(dynamic)
    for (int l = 0; l < NMAX; l++) {
        if (LMARK[l] == LMARKBULK) {
            for (int m = 0; m < LATTICE_VELOCITY_NUMBER; m++) {
                F[m][l] = FNEW[m][l];
            }
        }
    }
}

//Computes velocity field from F
void calcF2U(int l) {
    double fex = 0.0, fey = 0.0, density = 0.0;
  
    if (LMARK[l] == LMARKBULK) {
        for (int m = 0; m < LATTICE_VELOCITY_NUMBER; m++) {
            density += F[m][l];
            fex += F[m][l] * E[m][0];
            fey += F[m][l] * E[m][1];
        }

//        double factor = 0.5f;   //Correction due to discrete lattice effects
//        fex += factor * FORCE[1][l];
//        fey += factor * FORCE[2][l];
    
        U[1][l] = fex / density - MU * U[1][l] * DT / 2.0 / density;
        U[2][l] = fey / density - MU * U[2][l] * DT / 2.0 / density;
        U[0][l] = density;
    } else {
        U[0][l] = DENSITYINIT;
        U[1][l] = 0;
        U[2][l] = 0;
    }
}


//Compute stress tensor
void compute_sigma() {
    #pragma omp parallel for num_threads(STPROC) schedule(dynamic)
    for (int l = 0; l < NMAX; l++) {
        if (!(LMARK[l] == LMARKBULK)) continue;
        // NOTE: positive activity is extensile, negative activity is contractile
        // NOTE: Qxx = Q[0] and Qyy = -Q[0], Qxy = Qyx = Q[1], Hxx = H[0] and Hyy = -H[0], Hxy = Hyx = H[1]
        //xx component of the stress tensor
        SIGMA[0][l] = -LAMBDA * H[0][l] - ACTIVITY[l] * Q[0][l];
        //xy component of the stress tensor
        SIGMA[1][l] = -LAMBDA * H[1][l] + Q[0][l] * H[1][l] - Q[1][l] * H[0][l] - (H[0][l] * Q[1][l] - H[1][l] * Q[0][l]) - ACTIVITY[l] * Q[1][l];
        //yx component of the stress tensor
        SIGMA[2][l] = -LAMBDA * H[1][l] + Q[1][l] * H[0][l] - Q[0][l] * H[1][l] - (H[1][l] * Q[0][l] - H[0][l] * Q[1][l]) - ACTIVITY[l] * Q[1][l];
        //yy component of the stress tensor
        SIGMA[3][l] = +LAMBDA * H[0][l] + ACTIVITY[l] * Q[0][l];
    }

    //Impose open boundaries
    #pragma omp parallel for num_threads(STPROC) schedule(dynamic)
    for (int l = 0; l < NMAX; l++) {
        int lobc = calcLobc(l);
        if (lobc != l) {
            for (int m = 0; m < 2; m++) {
                SIGMA[m][l] = SIGMA[m][lobc];
            }
        }
    }
}
