//utility_functions.h
void initialiseE();
void write_velocity(int t);
void write_orientation(int t);

//lattice_Boltzmann.c
void compute_LB_step(double **U, int **E, double **FNEW, double **F, double **FEQ, double **P, double **Q, double **H, double **SIGMA, double *WKONST, double *ACTIVITY, char *LMARK);

__global__
void computeFeq(double **U, int **E, double **FEQ, double *WKONST);

__global__
void computeSigma(double **SIGMA, double **Q, double **H, double *ACTIVITY, char *LMARK);

__global__
void computeP(double **U, int **E, double **P, double **SIGMA, double *WKONST, char *LMARK);

__global__
void computeFNEW(double **FNEW, double **F, double **FEQ, int **E, double **P);

//__global__
//void enforceBoundaryConditionsAndCalcFNEW2F(double **FNEW, double **F, double **U);

//__global__
//void calcF2U(double **U, int **E, double **F, double **SIGMA, char *LMARK);

__global__
void enforceBCAndCalcFNEW2F2U(double **FNEW, double **F, double **U, int **E, double **SIGMA, char *LMARK);

//finite_difference.c
void compute_FD_step(double **U, double **Q, double **QNEW, double **H, char *LMARK);

__global__
void computeQNEW(double **U, double **Q, double **QNEW, double **H, char *LMARK);

__global__
void calcQNEW2Q(double **Q, double **QNEW, char *LMARK);
