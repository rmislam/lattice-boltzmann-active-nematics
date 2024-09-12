
//Make a step in the Q-tensor finite difference evolution
void compute_FD_step() {
    #pragma omp parallel for num_threads(STPROC) schedule(dynamic)
    for (int l = 0; l < NMAX; l++) {
        if (LMARK[l] == LMARKBULK) {
            double u1[2], u2[2], u3[2];
            //Elastic part
            double Q_laplacian[2];
            compute_Q_laplacian(l, Q_laplacian);
            compute_u1(l, Q_laplacian, u1);
            //Shear part
            compute_u2(l, u2);
            //Advective part
            compute_u3(l, u3);
            //Make the time step
            for (int m = 0; m < 2; m++) {
                QNEW[m][l] = Q[m][l] + DT * u1[m] + u2[m] + u3[m];
            }
        }
    }
    //Write QNEW back to Q
    calcQNEW2Q();
}

//Write QNEW back to Q and implement the periodic boundaries
void calcQNEW2Q() {
    #pragma omp parallel for num_threads(STPROC) schedule(dynamic)
    for (int l = 0; l < NMAX; l++) {
        int lpbc = calcLpbc(l);
        for(int m = 0; m < 2; m++) {
            Q[m][l] = QNEW[m][lpbc];
        }

        /*  // attempt at homeotropic anchoring condition
        if (j_vr(l) == 0 || j_vr(l) == J - 1) {
            double angle = 0.5  * M_PI;  // homeotropic anchoring condition at top and bottom walls (?)
            double degree_of_order = 1.;
            Q[0][l] = degree_of_order / 2.0 * cos(2 * angle);   //Qxx component
            Q[1][l] = degree_of_order / 2.0 * sin(2 * angle);   //Qxy component
        }
        */
    }
}

//Compute local Q tensor laplacian
void compute_Q_laplacian(int l, double* Q_laplacian) {
    int lxp = l + 1;    //Indices of neighbouring mesh points
    int lxm = l - 1;
	int lyp = l + I;
    int lym = l - I;
    
    for (int m = 0; m < 2; m++) {
        Q_laplacian[m] = Q[m][lxp] + Q[m][lxm] + Q[m][lyp] + Q[m][lym] - 4 * Q[m][l];
    }
}

//Compute the elastic contribution to the Q-tensor dynamics
void compute_u1(int l, double* Q_laplacian, double* u1) {
    double temp = 4 * (Q[0][l] * Q[0][l] + Q[1][l] * Q[1][l]);
    u1[0] = Q_laplacian[0] - 2 * C * Q[0][l] * (temp - 1.0);
    u1[1] = Q_laplacian[1] - 2 * C * Q[1][l] * (temp - 1.0);
    
    //Write u1 as H in a global matrix
    H[0][l] = u1[0];
    H[1][l] = u1[1];
}

//Compute the shear contribution to the Q-tensor dynamics
void compute_u2(int l, double* u2) {
    double uxx = (U[1][l + 1] - U[1][l - 1]) / 2.0;
    double uxy = 0.5 * ( (U[1][l + I] - U[1][l - I]) / 2.0 + (U[2][l + 1] - U[2][l - 1]) / 2.0 );
    u2[0] = LAMBDA * uxx;
    u2[1] = LAMBDA * uxy;
}

//Compute the advective contribution to the Q-tensor dynamics
void compute_u3(int l, double* u3) {
    double dQxxdx = (Q[0][l + 1] - Q[0][l - 1]) / 2.0;
    double dQxydx = (Q[1][l + 1] - Q[1][l - 1]) / 2.0;
    double dQxxdy = (Q[0][l + I] - Q[0][l - I]) / 2.0;
    double dQxydy = (Q[1][l + I] - Q[1][l - I]) / 2.0;
    
    u3[0] = -(U[1][l] * dQxxdx + U[2][l] * dQxxdy);
    u3[1] = -(U[1][l] * dQxydx + U[2][l] * dQxydy);
}
