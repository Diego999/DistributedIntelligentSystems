#include "zpso.h"

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <webots/robot.h>
#include <webots/emitter.h>
#include <webots/supervisor.h>

/* Types of fitness evaluations */
#define EVOLVE 0          // Find new fitness
#define EVOLVE_AVG 1      // Average new fitness into total
#define SELECT 2          // Find more accurate fitness for best selection

/* Size of swarm data must be global variables */
int swarmsize;
int datasize;
int robots;
int nb;
char label[20];

/* Particle swarm optimization function                                      */
/*                                                                           */
/* Parameters:                                                               */
/* n_swarmsize: number of particles in swarm                                 */
/* nb:          number of neighbors on each side of particle in neighborhood */
/* lweight:     max random value for local weight                            */
/* nbweight:    max random value for neighborhood weight                     */
/* vmax:        maximum velocity value                                       */
/* min:         minimum initial value of particle element                    */
/* max:         maximum initial value of particle element                    */
/* iterations:  number of iterations to run in the optimization              */
/* n_datasize:  number of elements in particle                               */
double* pso(int n_swarmsize, int n_nb, double lweight, double nbweight, double vmax, double min, double max, int iterations, int n_datasize, int n_robots) {
  double swarm[n_swarmsize][n_datasize];    // Swarm of particles
  double perf[n_swarmsize];                 // Current local performance of swarm
  double lbest[n_swarmsize][n_datasize];    // Current best local swarm
  double lbestperf[n_swarmsize];            // Current best local performance
  double lbestage[n_swarmsize];             // Life length of best local swarm
  double nbbest[n_swarmsize][n_datasize];   // Current best neighborhood
  double nbbestperf[n_swarmsize];           // Current best neighborhood performance
  double v[n_swarmsize][n_datasize];        // Preference indicator
  int neighbors[n_swarmsize][n_swarmsize];  // Neighbor matrix
  int i,j,k;                                // FOR-loop counters
  double bestperf;                          // Performance of evolved solution

    // Set global variables
    swarmsize = n_swarmsize;
    datasize = n_datasize;
    robots = n_robots;
    nb = n_nb;

    //printf("NOISY = %d\n",NOISY);
    sprintf(label, "Iteration: 0");
    wb_supervisor_set_label(0,label,0.01,0.01,0.1,0xffffff,0);
    // Seed the random generator
    srand(time(NULL));

    // Setup neighborhood
    for (i = 0; i < swarmsize; i++) {
        for (j = 0; j < swarmsize; j++) {
            if (mod(i-j+nb,swarmsize) <= 2*nb)
                neighbors[i][j] = 1;
            else
                neighbors[i][j] = 0;
        }
    }

    // Initialize the swarm
    for (i = 0; i < swarmsize; i++) {
        for (j = 0; j < datasize; j++) {
            // Randomly assign initial value in [min,max]
            swarm[i][j] = (max-min)*rnd()+min;
            lbest[i][j] = swarm[i][j];           // Best configurations are initially current configurations
            nbbest[i][j] = swarm[i][j];
            v[i][j] = 2.0*vmax*rnd()-vmax;         // Random initial velocity
        }
    }

    // Best performances are initially current performances
    printf("aa\n");
    findPerformance(swarm,perf,NULL,EVOLVE,robots,neighbors);
    printf("bb\n");
    for (i = 0; i < swarmsize; i++) {
        lbestperf[i] = perf[i];
        lbestage[i] = 1.0;                    // One performance so far
        nbbestperf[i] = perf[i];
    }
    updateNBPerf(lbest,lbestperf,nbbest,nbbestperf,neighbors);  // Find best neighborhood performances

    printf("Swarm initialized\n");

    // Run optimization
    for (k = 0; k < iterations; k++) {

    printf("Iteration %d\n",k);
    sprintf(label, "Iteration: %d",k+1);
    wb_supervisor_set_label(0,label,0.01,0.01,0.1,0xffffff,0);
    // Update preferences and generate new particles
    for (i = 0; i < swarmsize; i++) {
        for (j = 0; j < datasize; j++) {
            // Adjust preferences
            v[i][j] += lweight*rnd()*(lbest[i][j] - swarm[i][j]) + nbweight*rnd()*(nbbest[i][j] - swarm[i][j]);
            v[i][j] *= 0.6;
            swarm[i][j] += v[i][j];
        }
    }
    printf("a\n");
    // RE-EVALUATE PERFORMANCES OF PREVIOUS BESTS
    #if NOISY == 1
        findPerformance(lbest,lbestperf,lbestage,EVOLVE_AVG,robots,neighbors);
    #endif  
    printf("b\n");
    // Find new performance
    findPerformance(swarm,perf,NULL,EVOLVE,robots,neighbors);
    printf("c\n");
    // Update best local performance
    updateLocalPerf(swarm,perf,lbest,lbestperf,lbestage);
    printf("d\n");
    // Update best neighborhood performance
    updateNBPerf(lbest,lbestperf,nbbest,nbbestperf,neighbors);

    double temp[datasize];
    bestperf = bestResult(lbest,lbestperf,temp);
    printf("%f\n",bestperf);

    }

    // Find best result achieved
    double* best;
    best = malloc(sizeof(double)*datasize);
    findPerformance(lbest,lbestperf,NULL,SELECT,robots,neighbors);
    bestperf = bestResult(lbest,lbestperf,best);

    printf("Best performance found\n");
    printf("Performance: %f\n",bestperf);

    return best;
}

// Generate random number in [0,1]
double rnd(void) {
    return ((double)rand())/((double)RAND_MAX);
}

// Find the current performance of the swarm.
void findPerformance(double swarm[swarmsize][datasize], double perf[swarmsize], 
             double age[swarmsize], char type, int robots, 
             int neighbors[swarmsize][swarmsize]) {
    double particles[robots][datasize];
    double fit[robots];
    int i,j,k;                   // FOR-loop counters

    for (i = 0; i < swarmsize; i+=robots) {
        for (j=0;j<robots && i+j<swarmsize;j++) {
            for (k=0;k<datasize;k++) 
                particles[j][k] = swarm[i+j][k];
        }

        if (type == EVOLVE_AVG) {
            // Evalute current fitness
            fitness(particles,fit,neighbors);
            for (j=0;j<robots && i+j<swarmsize;j++) {
                perf[i+j] = (perf[i+j]*(age[i+j]-1)+fit[j])/age[i+j];
                age[i+j]++;
            }
        } else if (type == EVOLVE) {
            printf("ee\n");
        fitness(particles,fit,neighbors);
            for (j=0;j<robots && i+j<swarmsize;j++)
                perf[i+j] = fit[j];
        } else if (type == SELECT) {
            for (j=0;j<robots && i+j<swarmsize;j++)
                perf[i+j] = 0.0;
            for (k=0;k<5;k++) {
                fitness(particles,fit,neighbors);
                for (j=0;j<robots && i+j<swarmsize;j++)
                    perf[i+j] += fit[j];
            }
            for (j=0;j<robots && i+j<swarmsize;j++) {
                perf[i+j] /= 5.0;
            }
        }
    }
}

// Update the best performance of a single particle
void updateLocalPerf(double swarm[swarmsize][datasize], double perf[swarmsize], double lbest[swarmsize][datasize], double lbestperf[swarmsize], double lbestage[swarmsize]) {
  int i;                   // FOR-loop counters
  
    // If current performance of particle better than previous best, update previous best
    for (i = 0; i < swarmsize; i++) {
        if (perf[i] > lbestperf[i]) {
            copyParticle(lbest[i],swarm[i]);
            lbestperf[i] = perf[i];
            lbestage[i] = 1.0;
        }
    }
}

// Copy one particle to another
void copyParticle(double particle1[datasize], double particle2[datasize]) {
  int i;                   // FOR-loop counters

    // Copy one bit at a time
    for (i = 0; i < datasize; i++)
        particle1[i] = particle2[i];
}

// Update the best performance of a particle neighborhood
void updateNBPerf(double lbest[swarmsize][datasize], double lbestperf[swarmsize], 
              double nbbest[swarmsize][datasize], double nbbestperf[swarmsize], 
              int neighbors[swarmsize][swarmsize]) {
    int i,j;                   // FOR-loop counters

    // For each particle, check the best performances of its neighborhood (-NB to NB, with wraparound from swarmsize-1 to 0)
    for (i = 0; i < swarmsize; i++) {
        nbbestperf[i] = lbestperf[i];
        for (j = 0; j < swarmsize; j++) {
            // Make sure it's a valid particle
            if (!neighbors[i][j]) continue;
            // If current performance of particle better than previous best, update previous best
                if (lbestperf[j] > nbbestperf[i]) {
                    copyParticle(nbbest[i],lbest[j]);
                    nbbestperf[i] = lbestperf[j];
                }
        }
    }
}

// Find the modulus of an integer
int mod(int num, int base) {
    while (num >= base)
        num -= base;
    while (num < 0)      // Check for if number is negative to
        num += base;
    return num;
}


// S-function to transform v variable to [0,1]
double s(double v) {
    if (v > 5)
        return 1.0;
    else if (v < -5)
        return 0.0;
    else
        return 1.0/(1.0 + exp(-1*v));
}

// Find the best result found, set best to the particle, and return the performance
double bestResult(double lbest[swarmsize][datasize], double lbestperf[swarmsize], double best[datasize]) {
    double perf;         // Current best performance
    int i;               // FOR-loop counters

    // Start with the first particle as best
    copyParticle(best,lbest[0]);
    perf = lbestperf[0];

    // Iterate through the rest of the particles looking for better results
    for (i = 1; i < swarmsize; i++) {
        // If current performance of particle better than previous best, update previous best 
        if (lbestperf[i] > perf) {
            copyParticle(best,lbest[i]);
            perf = lbestperf[i];
        }
    }

    return perf;
}