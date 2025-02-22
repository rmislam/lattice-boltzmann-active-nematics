//Data file name
char FILE_NAME[80] = "output/active_nematic";

//Define M_PI
#define M_PI 3.14159265358979323846

//Number of mesh points
#define I 420
#define J 420

//Total number of points
#define NMAX (I*J)

//Total number of time steps in the simulation
#define TIME_STEPS 160000

//On how many steps the fields are written to a file
#define TIME_WRITE 2000

//On how many time steps to print progress
#define TIME_PRINT 5000

//How many time steps for Q relaxation
#define TIME_PRE_EVOL 10000

//How many FD steps to take for Q for each LB step
#define N_EVOL_Q 2

//Number of velocity vectors within the model -- in this case D2Q9 -- 9 velocities
#define LATTICE_VELOCITY_NUMBER 9

//Number of CPU cores for the paralelisation
int STPROC = 2;

//Free energy elastic constant
#define L 0.1 // 1.

//Constant A from free energy equation
#define A (0.1 * (1 - 3.5 / 3.0))

//Constant B from free energy equation
#define B (-1 * 0.1 * 3.5)

//Free energy phase parameter
#define C (0.1 * 3.5)
// 0.35 // 0.01 //0.1 // 1.  // -0.01

//Time step
#define DT 1 //0.05 // 0.002 // 0.0001   // try 0.01 or 0.005

//relaxation time of the Lb scheme
//physical kinematic viscosity = (1 / 3) * (TAUF - 0.5) * dx^2 / DT
#define TAUF 1 //(0.55 * DT) // corresponds to viscosity of 1/3, taking into account GAMMA // (1.0 * DT)  // must be greater than 0.5 * DT

//Density parameter
#define	DENSITYINIT 1 //(2.0 / DT)

//Molecular field coefficient
#define GAMMA 0.1

//Flow-aligning parameter
#define XI 0.8

//Ekman linear friction coefficient
#define MU 0.01

//Activity
#define ALPHA 0.0035 //0.022 // 1. //5.  // try 0.2 or 0.3 // 0.4 is working

//Logical markers for bulk and boundary points
#define LMARK_BULK 2          // channel bulk
#define LMARK_OBS_BULK 4      // obstacle bulk

#define LMARK_BOT_WALL 6
#define LMARK_TOP_WALL 8
#define LMARK_LEFT_WALL 10
#define LMARK_RIGHT_WALL 12

#define LMARK_CORNER_BOT_LEFT 14
#define LMARK_CORNER_BOT_RIGHT 16
#define LMARK_CORNER_TOP_LEFT 18
#define LMARK_CORNER_TOP_RIGHT 20

#define LMARK_BOT_OUTLET 22
#define LMARK_TOP_OUTLET 24
#define LMARK_LEFT_OUTLET 26
#define LMARK_RIGHT_OUTLET 28

//CUDA parameters
#define BLOCK_SIZE 32