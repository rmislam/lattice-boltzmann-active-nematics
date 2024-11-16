//Makes a main step in the LB part
void compute_LB_step() {
    computeFeq();   //Compute the equilibrium distribution
    computeP(); //Compute the forcing terms
    
    #pragma omp parallel for num_threads(STPROC) schedule(dynamic)
    for (int l = 0; l < NMAX; l++) {
        for (int m = 0; m < LATTICE_VELOCITY_NUMBER; m++) {
            int lnew = calcLBlnew(l, m);   //Streaming location
            if (lnew >= NMAX || lnew < 0) continue;
            FNEW[m][lnew] = F[m][l] + DT * ((FEQ[m][l] - F[m][l]) / TAUF + P[m][l]);
        }
    }
    
    enforceBoundaryConditions();  // NOTE: be careful about copying FNEW to F at boundaries
    calcFNEW2F();  //Copies FNEW to F
    
    //Compute the velocity field
    #pragma omp parallel for num_threads(STPROC) schedule(dynamic)
    for (int l = 0; l < NMAX; l++) calcF2U(l);
}

 //Computes equilibrium distribution function
void computeFeq() {
    #pragma omp parallel for num_threads(STPROC) schedule(dynamic)
    for (int l = 0; l < NMAX; l++) {
        double u2 = U[1][l] * U[1][l] + U[2][l] * U[2][l];	//Velocity squared

        for (int m = 0; m < LATTICE_VELOCITY_NUMBER; m++) {
            double ue = U[1][l] * E[m][0] + U[2][l] * E[m][1];		//Product of velocity and characteristic vectors
            FEQ[m][l] = WKONST[m] * U[0][l] * (1.0 + 3.0 * ue + 9.0 / 2.0 * ue * ue - 1.5 * u2);  // NOTE: no need to modify
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

        double uF = U[1][l] * forceX + U[2][l] * forceY;  //Product of force and velocity
        for (int m = 0; m < LATTICE_VELOCITY_NUMBER; m++) {
            double ue = U[1][l] * E[m][0] + U[2][l] * E[m][1];	//Product of velocity and characteristic vectors
            double eF = E[m][0] * forceX + E[m][1] * forceY;  //Product of force and characteristic vectors
            P[m][l] = (1.0 - DT / 2.0 / TAUF) * WKONST[m] * (3.0 * eF - 3.0 * uF + 9.0 * ue * eF);  // "source term" S_i from Kruger 6.2
        }
    }
}

//Computes streaming location
int calcLBlnew(int l, int m) {
    int i, j;
    i = i_vr(l) + E[m][0];
    j = j_vr(l) + E[m][1];
    return i + j * I;
}

//Enforce bounce back for no-slip velocity condition
void enforceBoundaryConditions() {
    #pragma omp parallel for num_threads(STPROC) schedule(dynamic)
    for (int l = 0; l < NMAX; l++) {
        if (j_vr(l) == 0 && !(LMARK[l] == LMARKOBSBULK || LMARK[l] == LMARKOBS_BOT_LEFT_CORNER || LMARK[l] == LMARKOBS_BOT_RIGHT_CORNER)) {  // bottom edge of channel -- no-slip
            FNEW[2][l] = F[4][l];
            FNEW[5][l] = F[7][l];
            FNEW[6][l] = F[8][l];
        } else if (j_vr(l) == 0 && LMARK[l] == LMARKOBSBULK) {
            FNEW[4][l] = F[2][l];
            FNEW[7][l] = F[5][l];
            FNEW[8][l] = F[6][l];
        } else if (LMARK[l] == LMARKOBSTOP) {  // top of obstacle -- no-slip
            FNEW[2][l] = F[4][l];
            FNEW[5][l] = F[7][l];
            FNEW[6][l] = F[8][l];
        } else if (j_vr(l) == J - 1) {  // top edge of channel -- no-slip
            FNEW[4][l] = F[2][l];
            FNEW[7][l] = F[5][l];
            FNEW[8][l] = F[6][l];
        } else if (LMARK[l] == LMARKOBSLEFT) {
            FNEW[6][l] = F[8][l];
            FNEW[3][l] = F[1][l];
            FNEW[7][l] = F[5][l];
        } else if (LMARK[l] == LMARKOBSRIGHT) {
            FNEW[5][l] = F[7][l];
            FNEW[1][l] = F[3][l];
            FNEW[8][l] = F[6][l];
        } else if (LMARK[l] == LMARKOBS_TOP_LEFT_CORNER) {
            FNEW[6][l] = F[8][l];
            FNEW[2][l] = F[4][l];
            FNEW[3][l] = F[1][l];
        } else if (LMARK[l] == LMARKOBS_TOP_RIGHT_CORNER) {
            FNEW[5][l] = F[7][l];
            FNEW[2][l] = F[4][l];
            FNEW[1][l] = F[3][l];
        } else if (LMARK[l] == LMARKOBS_BOT_LEFT_CORNER) {
            FNEW[6][l] = F[8][l];
        } else if (LMARK[l] == LMARKOBS_BOT_RIGHT_CORNER) {
            FNEW[5][l] = F[7][l];
        } else if (i_vr(l) == 0 && j_vr(l) > 0 && j_vr(l) < J - 1) {  // left edge -- fixed velocity inlet
            FNEW[1][l] = F[3][l] + 2.0 * U[0][l] * INLET_VELOCITY / 3.0;
            FNEW[5][l] = F[7][l] - 0.5 * (F[2][l] - F[4][l]) + U[0][l] * INLET_VELOCITY / 6.0;
            FNEW[8][l] = F[6][l] + 0.5 * (F[2][l] - F[4][l]) + U[0][l] * INLET_VELOCITY / 6.0;
        } else if (i_vr(l) == I - 1 && j_vr(l) > 0 && j_vr(l) < J - 1) {  // right edge -- open (absorbing) boundary
            //for (int m = 0; m < LATTICE_VELOCITY_NUMBER; m++) {
            //    FNEW[m][l] = F[m][l - 1];
            //}
            FNEW[3][l] = F[3][l - 1];
            FNEW[6][l] = F[6][l - 1];
            FNEW[7][l] = F[7][l - 1];
        }
    }
}

//Copies FNEW to F
void calcFNEW2F() {
    #pragma omp parallel for num_threads(STPROC) schedule(dynamic)
    for (int l = 0; l < NMAX; l++) {
        for (int m = 0; m < LATTICE_VELOCITY_NUMBER; m++) {
            F[m][l] = FNEW[m][l];
        }
    }
}

//Computes velocity field from F
void calcF2U(int l) {
    if (j_vr(l) == 0 || j_vr(l) == J - 1) {  // top and bottom edges -- no-slip condition
        U[1][l] = 0;
        U[2][l] = 0;
    } else if (i_vr(l) == 0 && j_vr(l) > 0 && j_vr(l) < J - 1) {  // left edge -- constant velocity (inlet) // NOTE: maybe ok to delete?
        U[1][l] = INLET_VELOCITY;
        U[2][l] = 0;
    } else if (LMARK[l] == LMARKOBSBULK || LMARK[l] == LMARKOBSLEFT || LMARK[l] == LMARKOBSRIGHT || LMARK[l] == LMARKOBSTOP || LMARK[l] == LMARKOBS_TOP_LEFT_CORNER || LMARK[l] == LMARKOBS_TOP_RIGHT_CORNER || LMARK[l] == LMARKOBS_BOT_LEFT_CORNER || LMARK[l] == LMARKOBS_BOT_RIGHT_CORNER) {
        U[1][l] = 0;
        U[2][l] = 0;
    //} else if (i_vr(l) == I - 1 && j_vr(l) > 0 && j_vr(l) < J - 1) {  // right edge -- open (absorbing) boundary //  NOTE: This ruins everything!
    //    U[0][l] = U[0][l - 1];
    //    U[1][l] = U[1][l - 1];
    //    U[2][l] = U[2][l - 1];
    } else if (LMARK[l] == LMARKBULK) {
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
    }
}

//Compute stress tensor
void compute_sigma() {
    #pragma omp parallel for num_threads(STPROC) schedule(dynamic)
    for (int l = 0; l < NMAX; l++) {
        if (!(LMARK[l] == LMARKBULK)) continue;
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
