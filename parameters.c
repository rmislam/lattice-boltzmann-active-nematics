//Data file name
char FILE_NAME[80] = "output/active_nematic";

//Define M_PI
#define M_PI 3.14159265358979323846

//Number of mesh points
#define I 121
#define J 121

//Total number of points
#define NMAX (I*J)

//Total number of time steps in the simulation
#define TIME_STEPS 16000 
//On how many steps the fields are written to a file
#define TIME_WRITE 50

//Number of velocity vectors within the model -- in this case D2Q9 -- 9 velocities
#define LATTICE_VELOCITY_NUMBER 9

//Number of CPU cores for the paralelisation
int STPROC = 2;

//Free energy phase parameter
#define C 1.0

//Time step
#define DT 0.01

//relaxation time of the Lb scheme
#define TAUF (2.0 * DT)

//Density parameter
#define	DENSITYINIT (2.0 / DT)

//Alignment parameter
#define LAMBDA 1.1  // 0.7

//Ekman linear friction coefficient
#define MU 0.01

//Activity
#define ALPHA 0.2

//Logical markers for bulk and boundary points
#define LMARKBULK 2    //bulk
#define LMARKBC 4      //boundary condition

//Activity pattern width and height fractions
#define AWIDTHFRAC 0.9
#define AHEIGHTFRAC 0.2
