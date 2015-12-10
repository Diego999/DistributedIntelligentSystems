#define NB_SENSOR   8
#define DATASIZE 2*(NB_SENSOR+2+1)         // Number of elements in particle
#define SWARMSIZE 10                    // Number of particles in swarm

// Functions
double* pso(int,int,double,double,double,double,double,int,int,int); // Run particle swarm optimization
void fitness(double[],double*,int[][SWARMSIZE]);                       // Fitness function for particle evolution
double rnd(void);                                                // Generate random number in [0,1]
void findPerformance(double[][DATASIZE],double[],double[],char,int,int[][SWARMSIZE]);  // Find the current performance of the swarm
void updateLocalPerf(double[][DATASIZE],double[],double[][DATASIZE],double[],double[]);   // Update the best performance of a single particle
void copyParticle(double[],double[]);                            // Copy value of one particle to another
void updateNBPerf(double[][DATASIZE],double[],double[][DATASIZE],double[],int[][SWARMSIZE]);  // Update the best performance of a particle neighborhood
int mod(int,int);                                                // Modulus function
double s(double);                                                // S-function to transform [-infinity,infinity] to [0,1]
double bestResult(double[][DATASIZE],double[],double[]);                 // Find the best result in a swarm
double rnd(void);