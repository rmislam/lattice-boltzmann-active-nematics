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
                QNEW[m][l] = Q[m][l] + DT * (GAMMA * u1[m] + u2[m] + u3[m]);
            }
        }

        double degree_of_order = 1.;
        double angle = 0.5 * M_PI; // vertical (homeotropic) anchoring at top and bottom boundaries  //0.25 * M_PI;

        // NOTE: maybe I can set this once at the beginning, and never set it again
        if (LMARK[l] == LMARKOBSLEFT || LMARK[l] == LMARKOBSRIGHT) {
            angle = 0; // horizontal (homeotropic) anchoring on left and right sides of obstacle
            QNEW[0][l] = degree_of_order / 2.0 * cos(2 * angle);   //Qxx component
            QNEW[1][l] = degree_of_order / 2.0 * sin(2 * angle);   //Qxy component
        } else if (LMARK[l] == LMARKOBS_TOP_LEFT_CORNER || LMARK[l] == LMARKOBS_BOT_LEFT_CORNER) {
            angle = 0.75 * M_PI;
            QNEW[0][l] = degree_of_order / 2.0 * cos(2 * angle);   //Qxx component
            QNEW[1][l] = degree_of_order / 2.0 * sin(2 * angle);   //Qxy component
        } else if (LMARK[l] == LMARKOBS_TOP_RIGHT_CORNER || LMARK[l] == LMARKOBS_BOT_RIGHT_CORNER) {
            angle = 0.25 * M_PI;
            QNEW[0][l] = degree_of_order / 2.0 * cos(2 * angle);   //Qxx component
            QNEW[1][l] = degree_of_order / 2.0 * sin(2 * angle);   //Qxy component
        } else if (j_vr(l) == 0 || j_vr(l) == J - 1 || LMARK[l] == LMARKOBSTOP) {
            // infinite homeotropic anchoring at top and bottom of domain, and top of obstacle
            QNEW[0][l] = degree_of_order / 2.0 * cos(2 * angle);   //Qxx component
            QNEW[1][l] = degree_of_order / 2.0 * sin(2 * angle);   //Qxy component
        }
        // ignore Q inside obstacle bulk
    }
    //Write QNEW back to Q
    calcQNEW2Q();
}

//Write QNEW back to Q and implement the open boundaries on left and right edges
void calcQNEW2Q() {
    #pragma omp parallel for num_threads(STPROC) schedule(dynamic)
    for (int l = 0; l < NMAX; l++) {
        if (i_vr(l) == 0) {
            for(int m = 0; m < 2; m++) {
                Q[m][l] = QNEW[m][l + 1];
                //Q[m][l] = QNEW[m][l + I - 2];
            }
        }
        else if (i_vr(l) == I - 1) {
            for(int m = 0; m < 2; m++) {
                Q[m][l] = QNEW[m][l - 1];
                //Q[m][l] = QNEW[m][l - (I - 2)];
            }
        } else {
            for(int m = 0; m < 2; m++) {
                Q[m][l] = QNEW[m][l];
            }
        }
    }
}

//Compute local Q tensor laplacian
void compute_Q_laplacian(int l, double* Q_laplacian) {
    int lxp = l + 1;    //Indices of neighbouring mesh points
    int lxm = l - 1;
	int lyp = l + I;
    int lym = l - I;
    
    // Laplacian of tensor is simply the Laplacian of each component
    for (int m = 0; m < 2; m++) {
        // five-point stencil finite-difference
        Q_laplacian[m] = Q[m][lxp] + Q[m][lxm] + Q[m][lyp] + Q[m][lym] - 4 * Q[m][l];
    }
}

//Compute the elastic contribution to the Q-tensor dynamics
void compute_u1(int l, double* Q_laplacian, double* u1) {
    double Qsq00 = Q[0][l] * Q[0][l] + Q[1][l] * Q[1][l];

    u1[0] = L * Q_laplacian[0] - A * Q[0][l] - 2 * C * Q[0][l] * Qsq00;
    u1[1] = L * Q_laplacian[1] - A * Q[1][l] - 2 * C * Q[1][l] * Qsq00;

    //Write u1 as H in a global matrix
    // Remember H is molecular field
    H[0][l] = u1[0];
    H[1][l] = u1[1];
}

//Compute the shear contribution to the Q-tensor dynamics
void compute_u2(int l, double* u2) {
    double uxx = (U[1][l + 1] - U[1][l - 1]) / 2.0;
    double uyy = -uxx;  // Remember that div u = 0 (incompressibility condition), so uxx = -uyy
    double uxy = (U[1][l + I] - U[1][l - I]) / 2.0;
    double uyx = (U[2][l + 1] - U[2][l - 1]) / 2.0;

    // u2 is S in the Beris-Edwards equation
    u2[0] = Q[1][l] * (uxy - uyx) + XI * (uxx + 2 * Q[0][l] * Q[0][l] * (uyy - uxx) - 2 * Q[0][l] * Q[1][l] * (uyx + uxy));
    u2[1] = (0.5 * XI - Q[0][l]) * uxy + (0.5 * XI + Q[0][l]) * uyx;
}

//Compute the advective contribution to the Q-tensor dynamics
void compute_u3(int l, double* u3) {
    // Gradient of Q
    double dQxxdx = (Q[0][l + 1] - Q[0][l - 1]) / 2.0;
    double dQxydx = (Q[1][l + 1] - Q[1][l - 1]) / 2.0;
    double dQxxdy = (Q[0][l + I] - Q[0][l - I]) / 2.0;
    double dQxydy = (Q[1][l + I] - Q[1][l - I]) / 2.0;
    
    // -u dot grad Q
    u3[0] = -(U[1][l] * dQxxdx + U[2][l] * dQxxdy);
    u3[1] = -(U[1][l] * dQxydx + U[2][l] * dQxydy);
}
