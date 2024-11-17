//utility_functions.h
int i_vr(int l);
int j_vr(int l);
void initialiseE();
void write_velocity(int t);
void write_orientation(int t);

//lattice_Boltzmann.c
void computeFeq();
void compute_LB_step();
void computeP();
void computeBounceBack();
int calcLBlnew(int l, int m);
void calcFNEW2F();
void calcF2U(int l);
void compute_sigma();
void enforceBoundaryConditions();

//finite_difference.c
void compute_FD_step();
void calcQNEW2Q();
void compute_Q_laplacian(int l, double* Q_laplacian);
void compute_u1(int l,double* Q_laplacian,double* u1);
void compute_u2(int l, double* u2);
void compute_u3(int l, double* u3);

