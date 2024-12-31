//Make a step in the Q-tensor finite difference evolution
void compute_FD_step(double **U, double **Q, double **QNEW, double **H, char *LMARK) {
    int num_blocks = (NMAX + BLOCK_SIZE - 1) / BLOCK_SIZE;

    computeQNEW<<<num_blocks, BLOCK_SIZE>>>(U, Q, QNEW, H, LMARK);
    cudaDeviceSynchronize();

    calcQNEW2Q<<<num_blocks, BLOCK_SIZE>>>(Q, QNEW);
    cudaDeviceSynchronize();
}

__global__
void computeQNEW(double **U, double **Q, double **QNEW, double **H, char *LMARK) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;

    for (int l = index; l < NMAX; l += stride) {
        if (LMARK[l] == LMARKBULK) {
            double u1[2], u2[2], u3[2];
            //Elastic part
            double Q_laplacian[2];


            // compute_Q_laplacian
            int lxp = l + 1;    //Indices of neighbouring mesh points
            int lxm = l - 1;
            int lyp = l + I;
            int lym = l - I;
            
            // Laplacian of tensor is simply the Laplacian of each component
            for (int m = 0; m < 2; m++) {
                // five-point stencil finite-difference
                Q_laplacian[m] = Q[m][lxp] + Q[m][lxm] + Q[m][lyp] + Q[m][lym] - 4 * Q[m][l];
            }


            // compute_u1
            double Qsq00 = Q[0][l] * Q[0][l] + Q[1][l] * Q[1][l];

            u1[0] = L * Q_laplacian[0] - A * Q[0][l] - 2 * C * Q[0][l] * Qsq00;
            u1[1] = L * Q_laplacian[1] - A * Q[1][l] - 2 * C * Q[1][l] * Qsq00;

            //Write u1 as H in a global matrix
            // Remember H is molecular field
            H[0][l] = u1[0];
            H[1][l] = u1[1];


            // compute_u2 (Shear part)
            double uxx = (U[1][l + 1] - U[1][l - 1]) / 2.0;
            double uyy = -uxx;  // Remember that div u = 0 (incompressibility condition), so uxx = -uyy
            double uxy = (U[1][l + I] - U[1][l - I]) / 2.0;
            double uyx = (U[2][l + 1] - U[2][l - 1]) / 2.0;

            // u2 is S in the Beris-Edwards equation
            u2[0] = Q[1][l] * (uxy - uyx) + XI * (uxx + 2 * Q[0][l] * Q[0][l] * (uyy - uxx) - 2 * Q[0][l] * Q[1][l] * (uyx + uxy));
            u2[1] = (0.5 * XI - Q[0][l]) * uxy + (0.5 * XI + Q[0][l]) * uyx;


            // compute_u3 (Advective part)
            // Gradient of Q
            double dQxxdx = (Q[0][l + 1] - Q[0][l - 1]) / 2.0;
            double dQxydx = (Q[1][l + 1] - Q[1][l - 1]) / 2.0;
            double dQxxdy = (Q[0][l + I] - Q[0][l - I]) / 2.0;
            double dQxydy = (Q[1][l + I] - Q[1][l - I]) / 2.0;
            
            // -u dot grad Q
            u3[0] = -(U[1][l] * dQxxdx + U[2][l] * dQxxdy);
            u3[1] = -(U[1][l] * dQxydx + U[2][l] * dQxydy);


            //Make the time step
            for (int m = 0; m < 2; m++) {
                QNEW[m][l] = Q[m][l] + DT * (GAMMA * u1[m] + u2[m] + u3[m]);
            }
        }

        double degree_of_order = 1.;
        double angle = 0.5 * M_PI; // vertical (homeotropic) anchoring at top and bottom boundaries  //0.25 * M_PI;

        // NOTE: maybe I can set this once at the beginning, and never set it again
        if (((l % (I * J)) / I) == 0 || ((l % (I * J)) / I) == J - 1) {
            // infinite homeotropic anchoring at top and bottom
            QNEW[0][l] = degree_of_order / 2.0 * cos(2 * angle);   //Qxx component
            QNEW[1][l] = degree_of_order / 2.0 * sin(2 * angle);   //Qxy component
        }
    }
}

//Write QNEW back to Q and implement the open boundaries on left and right edges
__global__
void calcQNEW2Q(double **Q, double **QNEW) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;

    for (int l = index; l < NMAX; l += stride) {
        if ((l % I) == 0) {
            for(int m = 0; m < 2; m++) {
                Q[m][l] = QNEW[m][l + 1];
                //Q[m][l] = QNEW[m][l + I - 2];
            }
        }
        else if ((l % I) == I - 1) {
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
