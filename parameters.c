//Data file name
char FILE_NAME[80] = "output/active_nematic";

//Define M_PI
#define M_PI 3.14159265358979323846

//Number of mesh points
#define I 101 //151 // 301
#define J 101 //41  // 81

//Total number of points
#define NMAX (I*J)

//Total number of time steps in the simulation
#define TIME_STEPS 300000
//On how many steps the fields are written to a file
#define TIME_WRITE 1000

// How many steps to wait before writing data files
#define TIME_START_WRITE 40000

// How long to wait for thermal equilibrium before applying defects and activity
#define WARM_UP_STEPS 40000

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
#define DT 0.035 // 0.002 // 0.0001   // try 0.01 or 0.005

//relaxation time of the Lb scheme
#define TAUF (0.55 * DT) // corresponds to viscosity of 1/3, taking into account GAMMA // (1.0 * DT)  // must be greater than 0.5 * DT

//Density parameter
#define	DENSITYINIT (2.0 / DT)

//Inlet velocity from left edge
#define INLET_VELOCITY 0.005

//Molecular field coefficient
#define GAMMA 0.1

//Flow-aligning parameter
#define XI 0.7

//Ekman linear friction coefficient
#define MU 0.01

//Activity alpha = ALPHA * L
#define ALPHA 0.22 // 1. //5.  // try 0.2 or 0.3 // 0.4 is working

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

//Activity pattern width and height fractions
#define AWIDTHFRAC 0.9
#define AHEIGHTFRAC 0.2
