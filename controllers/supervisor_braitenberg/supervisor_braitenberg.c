#include "zpso.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <webots/robot.h>
#include <webots/emitter.h>
#include <webots/receiver.h>
#include <webots/supervisor.h>


#define FLOCK_SIZE				5 		// Number of robots in flock
#define TIME_STEP				64		// [ms] Length of time step

WbNodeRef robs[FLOCK_SIZE];		// Robots nodes
WbFieldRef robs_trans[FLOCK_SIZE];	// Robots translation fields
WbFieldRef robs_rotation[FLOCK_SIZE];	// Robots rotation fields

float loc[FLOCK_SIZE][3];
float loc_old[FLOCK_SIZE][3];

#define RULE1_THRESHOLD 		0.2
#define fit_cluster_ref 		0.03
#define fit_orient_ref 			1.0
#define MAXSPEED                0.1287
#define SIZE_MAX_CLUSTER        FLOCK_SIZE+1
#define NB_MAX_CLUSTER          FLOCK_SIZE
#define END_CLUSTER_IDX         SIZE_MAX_CLUSTER
#define DH                      0.01
#define H_MAX                   5.00
//weights for the performance calculation
#define W_O    1.0   //orientation
#define W_C    1.0   //cohesion
#define W_V    1.0   //velocity
#define W_S    1.0   //entropy

int clusters[NB_MAX_CLUSTER][SIZE_MAX_CLUSTER]; //+1 because we want to add END_CLUSTER_IDX
float dist_robot[FLOCK_SIZE][FLOCK_SIZE];

int offset;					// Offset of robots number
float migrx, migrz;			// Migration vector
float orient_migr; 			// Migration orientation
int t;

/* PSO definitions */
#define NB 1                            // Number of neighbors on each side
#define LWEIGHT 2.0                     // Weight of attraction to personal best
#define NBWEIGHT 2.0                    // Weight of attraction to neighborhood best
#define VMAX 20.0                       // Maximum velocity particle can attain
#define MININIT -20.0                   // Lower bound on initialization value
#define MAXINIT 20.0                    // Upper bound on initialization value
#define ITS 10                          // Number of iterations to run
#define MAX_ROB FLOCK_SIZE
#define ROBOTS FLOCK_SIZE

#define NB_SENSOR   8
#define DATASIZE 2*NB_SENSOR+6          // Number of elements in particle
#define SWARMSIZE 30                    // Number of particles in swarm

/* Neighborhood types */
#define STANDARD    -1
#define RAND_NB      0
#define NCLOSE_NB    1
#define FIXEDRAD_NB  2

/* Fitness definitions */
#define FIT_ITS 180                     // Number of fitness steps to run during evolution

#define FINALRUNS 10
#define NEIGHBORHOOD STANDARD
#define RADIUS 0.8

WbDeviceTag emitter[MAX_ROB];
WbDeviceTag rec[MAX_ROB];
const double *loc_[MAX_ROB];
const double *rot[MAX_ROB];
double new_loc[MAX_ROB][3];
double new_rot[MAX_ROB][4];

void calc_fitness(double[DATASIZE],double*,int,int);
void random_pos(int);
void nRandom(int[][SWARMSIZE],int);
void nClosest(int[][SWARMSIZE],int);
void fixedRadius(int[][SWARMSIZE],double);

/*
 * Initialize flock position and devices
 */
void reset(void) {
	wb_robot_init();

	char rob[7] = "epuck0";
    char em[] = "emitter0";
    char receive[] = "receiver0";

	int i;
	for (i=0;i<FLOCK_SIZE;i++) {
        sprintf(rob,"epuck%d",i+offset);
        robs[i] = wb_supervisor_node_get_from_def(rob);
        robs_trans[i] = wb_supervisor_node_get_field(robs[i],"translation");
        loc_[i] = wb_supervisor_field_get_sf_vec3f(robs_trans[i]);
        new_loc[i][0] = loc_[i][0]; new_loc[i][1] = loc_[i][1]; new_loc[i][2] = loc_[i][2];
        robs_rotation[i] = wb_supervisor_node_get_field(robs[i],"rotation");
        rot[i] = wb_supervisor_field_get_sf_rotation(robs_rotation[i]);
        new_rot[i][0] = rot[i][0]; new_rot[i][1] = rot[i][1]; new_rot[i][2] = rot[i][2]; new_rot[i][3] = rot[i][3];
        emitter[i] = wb_robot_get_device(em);
        if (emitter[i]==0) 
            printf("missing emitter %d\n",i);
        rec[i] = wb_robot_get_device(receive);
        em[7]++;
        receive[8]++;
    }
}

/*
 * Compute performance metric.
 */
void compute_fitness(float* fit_c, float* fit_o) {
	*fit_c = 0; *fit_o = 0;
	// Compute performance indices
	// Based on distance of the robots compared to the threshold and the deviation from the perfect angle towards
	// the migration goal
	float angle_diff;
	int i; int j;
	for (i=0;i<FLOCK_SIZE;i++) {
		for (j=i+1;j<FLOCK_SIZE;j++) 
			// Distance measure for each pair ob robots
			*fit_c += fabs(sqrtf(powf(loc[i][0]-loc[j][0],2) + powf(loc[i][1]-loc[j][1],2))-RULE1_THRESHOLD*2);
		
		// Angle measure for each robot
		angle_diff = fabsf(loc[i][2] - orient_migr);
		*fit_o += angle_diff > M_PI ? 2.0*M_PI-angle_diff : angle_diff;
	}
	*fit_c /= FLOCK_SIZE*(FLOCK_SIZE+1.0)/2.0;
	*fit_o /= ((float)FLOCK_SIZE);
}

/*
 * Compute performance metric orientation
 */
void compute_fitness_O(float* fit) {
    
    float angle_diff;
    float sum_cos = 0.0;
    float sum_sin = 0.0;
    
    int i;
    for (i=0; i < FLOCK_SIZE; i++){

        // compute angle with respect to migration angle
        angle_diff = fabsf(loc[i][2] - orient_migr);
        angle_diff = angle_diff > M_PI ? 2*M_PI-angle_diff : angle_diff;
        
        // compute cos of the angle
        sum_cos += cos(angle_diff);
        
        // compute sin of the angle
        sum_sin += sin(angle_diff);
    }
    
    // compute the orientation metric
    *fit = (1/ (double) FLOCK_SIZE) * sqrtf( powf(sum_cos, 2) + powf(sum_sin, 2));
}

/*
 * Compute performance metric orientation
 */
void compute_fitness_C(float* fit) {
    
    float x_bar = 0.0;
    float z_bar = 0.0;
    float dist  = 0.0;
    
    int i;
    
    // compute center of mass
    for (i=0; i < FLOCK_SIZE; i++){
        x_bar += (1/ (double) FLOCK_SIZE) * loc[i][0];
        z_bar += (1/ (double) FLOCK_SIZE) * loc[i][1];
    }
    
    // compute the sum distance
    for (i=0; i < FLOCK_SIZE; i++){
        dist += sqrtf( powf(loc[i][0] - x_bar, 2) + powf(loc[i][1] - z_bar, 2));
    }
    
    // compute the metric
    *fit = 1/ (1 + (1/ (double) FLOCK_SIZE)*dist);
    
}

/*
 * Compute performance metric orientation
 */
void compute_fitness_V(float* fit) {
    
    int i;
    
    // initialize current center of max coordinate
    float x_bar = 0.0;
    float z_bar = 0.0;
    
    // initialize previous center of max coordinate
    float x_bar_old = 0.0;
    float z_bar_old = 0.0;
    
    // initialize migration urge direction components
    float psi_x;
    float psi_z;
    
    // initialize scalar projection of dx onto psi
    float proj;
    
    for (i=0; i < FLOCK_SIZE; i++){
        // compute the current center of mass
        x_bar += (1/ (double) FLOCK_SIZE) * loc[i][0];
        z_bar += (1/ (double) FLOCK_SIZE) * loc[i][2];
        // compute the old previous center of mass
        x_bar_old += (1/ (double) FLOCK_SIZE) * loc_old[i][0];
        z_bar_old += (1/ (double) FLOCK_SIZE) * loc_old[i][2];
    }
    
    // compute migration urge components
    psi_x = cos(orient_migr);
    psi_z = sin(orient_migr);
    
    // compute the projection
    proj = (x_bar - x_bar_old)*(psi_x) + (z_bar - z_bar_old)*psi_z;
    
    //printf("%f\n", x_bar);
    //printf("%f\n", z_bar);
    //printf("%f\n", x_bar_old);
    //printf("%f\n", z_bar_old);
    //printf("%f\n", proj);
    
    // compute the metric
    *fit = (1.0 / MAXSPEED) * proj;
    
}

int compare (const void * a, const void * b)
{
    return ( *(int*)a - *(int*)b );
}

float compute_H(float h) {
    float r = h/2.0;
    int i,j,k;
    int temp;
    int cluster_idx_robot = 0;
    int cluster_idx = 0;
    float x, y;
    float res = 0;
    float pk = 0;
    int size = 0;

    // Initialize cluster array
    for(i = 0; i < NB_MAX_CLUSTER; ++i)
        for(j = 0; j < SIZE_MAX_CLUSTER; ++j)
            clusters[i][j] = END_CLUSTER_IDX;

    // Compute distance matrix
    for(i = 0; i < FLOCK_SIZE; ++i)
        for(j = 0; j < FLOCK_SIZE; ++j) {
            x = loc[i][0]-loc[j][0];
            y = loc[i][1]-loc[j][1];

            dist_robot[i][j] = sqrtf(x*x + y*y);
            }

    // Find clusters for every robots
    for(i = 0; i < FLOCK_SIZE; ++i) {
        // Every robot is a cluster at least of itself
        cluster_idx_robot = 0;
        clusters[cluster_idx][cluster_idx_robot] = i;

        // Find other robots to add to the cluster
        for(j = 0; j < FLOCK_SIZE; ++j)
            if(i != j && dist_robot[i][j] <= r)
                clusters[cluster_idx][++cluster_idx_robot] = j;

        ++cluster_idx;
    }

    // Sort all clusters
    for(i = 0; i < NB_MAX_CLUSTER && clusters[i][0] != END_CLUSTER_IDX; ++i)
        qsort(clusters[i], SIZE_MAX_CLUSTER, sizeof(int), compare);

    // Keep only unique clusters
    for(i = 0; i < NB_MAX_CLUSTER; ++i) {
        // Avoid empty clusters
        if(clusters[i][0] == END_CLUSTER_IDX)
            continue;

        for(j = i+1; j < NB_MAX_CLUSTER; ++j) {
            // Avoid empty clusters
            if(clusters[j][0] == END_CLUSTER_IDX)
                continue;

            temp = 1;
            for(k = 0; temp == 1 && k < SIZE_MAX_CLUSTER; ++k)
                if(clusters[i][k] != clusters[j][k])
                    temp = 0;

            // We have find the same clusters, we remove it
            if(temp == 1)
                clusters[j][0] = END_CLUSTER_IDX;
        }
    }

    // Compute entropy
    res = 0.0f;
    for(i = 0; i < NB_MAX_CLUSTER; ++i)
        if(clusters[i][0] != END_CLUSTER_IDX) {
            size = 0;
            for(j = 0; j < SIZE_MAX_CLUSTER && clusters[i][j] != END_CLUSTER_IDX; ++j)
                ++size;
            
            pk = ((float)size)/((float)FLOCK_SIZE);
            res += (-pk * log2(pk));
        }

    return res;
}

void compute_fitness_S(float* fit) {
    float res = 0.0f;
    float h = DH;
    const float max_entropy = 500.0f;

    for(h = DH; h < H_MAX; h += DH)
        res += compute_H(h);
    
    //have an increasing entropy of [0, 1] (1 is best)
    res = ( max_entropy-res < 0)? 0.0f : ( max_entropy-res )/max_entropy ;
    
    *fit = res;
}

float instant_perf(){
   float fit_O = 0.0f;
   float fit_C = 0.0f;
   float fit_V = 0.0f;
   float fit_S = 0.0f;
   
   // compute the orientation metric
   compute_fitness_O(& fit_O);
   
   // compute the cohesion metric
   compute_fitness_C(& fit_C);
   
   // compute the velocity metric
   compute_fitness_V(& fit_V);

   // compute entropy metric
   compute_fitness_S(& fit_S);

   //printf("Result %f * %f * %f * %f = %f\n", W_O*fit_O , W_C*fit_C , W_V*fit_V , W_S*fit_S, W_O*fit_O * W_C*fit_C * W_V*fit_V * W_S*fit_S);
   return W_O*fit_O * W_C*fit_C * W_V*fit_V * W_S*fit_S;
}

/*
 * Main function.
 */
 
int main(int argc, char *args[]) {
    double *weights;                         // Evolved result
    double buffer[255];
    int i,j,k;

    wb_robot_init();
    printf("Particle Swarm Optimization Super Controller\n");
    reset();
    wb_robot_step(256);

    double fit, endfit, w[DATASIZE], f;
    double bestfit, bestw[DATASIZE];

    for(i = 0; i < FLOCK_SIZE; ++i)
    	wb_receiver_enable(rec[i],32);

    // Evolve controllers
    endfit = 0.0;
    bestfit = 0.0;
    for (j=0;j<10;j++) {
	    
        /* Get result of evolution */
        weights = pso(SWARMSIZE,NB,LWEIGHT,NBWEIGHT,VMAX,MININIT,MAXINIT,ITS,DATASIZE,ROBOTS);
        
        /* Calculate performance */
        fit = 0.0;
        for (k=0;k<DATASIZE;k++)
            w[k] = weights[k];
        
        calc_fitness(w,&f,FIT_ITS,MAX_ROB);
        fit = f;

        /* Check for new best fitness */
        if (fit > bestfit) {
            bestfit = fit;
        for (i = 0; i < DATASIZE; i++)
            bestw[i] = weights[i];
        }
        
        //printf("Performance: %.3f\n",fit);
        endfit += fit/10;

	    /* Send best controller to robots */
	    for (i=0;i<DATASIZE;i++) {
	        buffer[i] = bestw[i];
	    }
    }
    //printf("Average performance: %.3f\n",endfit);

    return 0;
}


// Randomly position specified robot
void random_pos(int rob_id) {
	wb_supervisor_field_set_sf_vec3f(wb_supervisor_node_get_field(robs[rob_id],"translation"), new_loc[rob_id]);
	wb_supervisor_field_set_sf_rotation(wb_supervisor_node_get_field(robs[rob_id],"rotation"), new_rot[rob_id]);
}

// Distribute fitness functions among robots
void calc_fitness(double weights[DATASIZE], double* fit, int its, int numRobs) {
    double buffer[255];
    int i;
    double perf = 0;
    int nbMeasure = 0;
    double* rbuffer;
    double local_perf = 0;

    /* Send data to robots */
    for(i = 0; i < FLOCK_SIZE; ++i)
	   	random_pos(i);

	/* Send best controller to robots */
	for (i=0;i<DATASIZE;i++) {
	    buffer[i] = weights[i];
	}
	buffer[DATASIZE] = 10000;
    for (i=0;i<ROBOTS;i++) {
        wb_emitter_send(emitter[i],(void *)buffer,(DATASIZE+1)*sizeof(double));
    }


	  /* Wait for response */
	  while (wb_receiver_get_queue_length(rec[0]) == 0) {
	  	for (i=0;i<FLOCK_SIZE;i++) {
	        // initialize old position
	        loc_old[i][0] = loc[i][0];
	        loc_old[i][1] = loc[i][1];
	        loc_old[i][2] = loc[i][2];
	        // initialize current position
		    loc[i][0] = wb_supervisor_field_get_sf_vec3f(robs_trans[i])[0];
		    loc[i][1] = wb_supervisor_field_get_sf_vec3f(robs_trans[i])[2];
		    loc[i][2] = wb_supervisor_field_get_sf_rotation(robs_rotation[i])[3];
	    }
	    local_perf = instant_perf();
	    if(local_perf > 0) {
		  	perf += local_perf;
		  	++nbMeasure;
		  }
	    wb_robot_step(64);
	}

	  /* Get fitness values */
	  for (i=0;i<FLOCK_SIZE;i++) {
	    rbuffer = (double *)wb_receiver_get_data(rec[i]);
	    wb_receiver_next_packet(rec[i]);
	  }

    *fit = perf / nbMeasure;
}

    // Evolution fitness function
void fitness(double weights[], double* fit, int neighbors[SWARMSIZE][SWARMSIZE]) {
    
    calc_fitness(weights,fit,FIT_ITS,ROBOTS);
    #if NEIGHBORHOOD == RAND_NB
        nRandom(neighbors,2*NB);
    #endif
    #if NEIGHBORHOOD == NCLOSE_NB
        nClosest(neighbors,2*NB);
    #endif
    #if NEIGHBORHOOD == FIXEDRAD_NB
        fixedRadius(neighbors,RADIUS);
    #endif
}

/* Get distance between robots */
double robdist(int i, int j) {
    return sqrt(pow(loc[i][0]-loc[j][0],2) + pow(loc[i][2]-loc[j][2],2));
}

/* Choose n random neighbors */
void nRandom(int neighbors[SWARMSIZE][SWARMSIZE], int numNB) {
    int i,j;

    /* Get neighbors for each robot */
    for (i = 0; i < ROBOTS; i++) {

    /* Clear old neighbors */
    for (j = 0; j < ROBOTS; j++)
        neighbors[i][j] = 0;

    /* Set new neighbors randomly */
    for (j = 0; j < numNB; j++)
        neighbors[i][(int)(SWARMSIZE*rnd())] = 1;
    }
}

/* Choose the n closest robots */
void nClosest(int neighbors[SWARMSIZE][SWARMSIZE], int numNB) {

    int r[numNB];
    int tempRob;
    double dist[numNB];
    double tempDist;
    int i,j,k;

    /* Get neighbors for each robot */
    for (i = 0; i < ROBOTS; i++) {

        /* Clear neighbors */
        for (j = 0; j < numNB; j++)
            dist[j] = 100000000;

        /* Find closest robots */
        for (j = 0; j < ROBOTS; j++) {

            /* Don't use self */
            if (i == j)
                continue;

            /* Check if smaller distance */
            if (dist[numNB-1] > robdist(i,j)) {
                dist[numNB-1] = robdist(i,j);
                r[numNB-1] = j;

                /* Move new distance to proper place */
                for (k = numNB-1; k > 0 && dist[k-1] > dist[k]; k--) {
                    tempDist = dist[k];
                    dist[k] = dist[k-1];
                    dist[k-1] = tempDist;
                    tempRob = r[k];
                    r[k] = r[k-1];
                    r[k-1] = tempRob;
                }
            }
        }

        /* Update neighbor table */
        for (j = 0; j < ROBOTS; j++)
            neighbors[i][j] = 0;
        for (j = 0; j < numNB; j++)
            neighbors[i][r[j]] = 1;
    }

}

/* Choose all robots within some range */
void fixedRadius(int neighbors[SWARMSIZE][SWARMSIZE], double radius) {
    int i,j;

    /* Get neighbors for each robot */
    for (i = 0; i < ROBOTS; i++) {
        /* Find robots within range */
        for (j = 0; j < ROBOTS; j++) {
            if (i == j) 
                continue;
            if (robdist(i,j) < radius) 
                neighbors[i][j] = 1;
            else 
                neighbors[i][j] = 0;
        }
    }
}

void step_rob() {
    wb_robot_step(64);
}
