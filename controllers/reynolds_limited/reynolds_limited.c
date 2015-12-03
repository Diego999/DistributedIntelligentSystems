#include <stdio.h>
#include <math.h>
#include <string.h>

#include <webots/robot.h>

#include <webots/differential_wheels.h>
#include <webots/distance_sensor.h>
#include <webots/emitter.h>
#include <webots/receiver.h>

//Braitenberg parameters for obstacle avoidance
#define NB_SENSORS      8    // Number of distance sensors
#define MIN_SENS        350       // Minimum sensibility value
#define MAX_SENS        4096      // Maximum sensibility value
#define MAX_SPEED       800       // Maximum speed
#define FLOCK_SIZE      5    // Size of flock
#define TIME_STEP       64   // [ms] Length of time step

#define AXLE_LENGTH     0.052   // Distance between wheels of robot (meters)
#define SPEED_UNIT_RADS 0.00628    // Conversion factor from speed unit to radian per second
#define WHEEL_RADIUS    0.0205     // Wheel radius (meters)
#define DELTA_T         0.064   // Timestep (seconds)

#define MAX_TIMEOUT     1000 //# of timestamp before we forget the position of a robot (we kick it out of the flock)

#define RULE1_THRESHOLD   0.2    // Threshold to activate aggregation rule. default 0.20
#define RULE1_WEIGHT      0.3  // Weight of aggregation rule. default 0.20

#define RULE2_THRESHOLD   0.1    // Threshold to activate dispersion rule. default 0.1
#define RULE2_WEIGHT      1.0  // Weight of dispersion rule. default 1.0

#define RULE3_WEIGHT      0.01   // Weight of consistency rule. default 0.01

#define MIGRATION_WEIGHT  0.01   // Wheight of attraction towards the common goal. default 0.01

//weigths for the braientberg
int braiten_weight[16] = {17,  29,  34,  10, 8,  -38,-56, -76,
                          -72, -58, -36, 8,  10, 36,  28, 18 };

WbDeviceTag ds[NB_SENSORS];   // Handle for the infrared distance sensors
WbDeviceTag receiver2;     // Handle for the receiver node
WbDeviceTag emitter2;      // Handle for the emitter node

int robot_id_u, robot_id;  // Unique and normalized (between 0 and FLOCK_SIZE-1) robot ID

float relative_pos[FLOCK_SIZE][3];  // relative X, Z, Theta of all robots
float prev_relative_pos[FLOCK_SIZE][3];   // Previous relative  X, Z, Theta values
float my_position[3];         // X, Z, Theta of the current robot
float prev_my_position[3];       // X, Z, Theta of the current robot in the previous time step
float speed[FLOCK_SIZE][2];      // Speeds calculated with Reynold's rules
float relative_speed[FLOCK_SIZE][2];   // Speeds calculated with Reynold's rules
int initialized[FLOCK_SIZE];     // != 0 if initial positions have been received
float migr[2] = {0,-50};           // Migration vector
char* robot_name;
char robot_number[8];

unsigned int timestamp[FLOCK_SIZE];

int maxTimestamp;

/*
 * Reset the robot's devices and get its ID
 */
static void reset() 
{
   int i;
   wb_robot_init();

   receiver2 = wb_robot_get_device("receiver2");
   emitter2 = wb_robot_get_device("emitter2");
   
   char s[4]="ps0";
   for(i=0; i<NB_SENSORS;i++) {
      ds[i]=wb_robot_get_device(s); // the device name is specified in the world file
      s[2]++;           // increases the device number
   }
   
   for(i=0; i<FLOCK_SIZE; i++) {
      timestamp[i] = 0;
   }
   maxTimestamp = 1;
   
   robot_name=(char*) wb_robot_get_name();
   
   for(i=0;i<NB_SENSORS;i++)
      wb_distance_sensor_enable(ds[i],64);

   wb_receiver_enable(receiver2,64);

   //Reading the robot's name. Pay attention to name specification when adding robots to the simulation!
   sscanf(robot_name,"epuck%d",&robot_id_u); // read robot id from the robot's name
   robot_id = robot_id_u%FLOCK_SIZE;     // normalize between 0 and FLOCK_SIZE-1
   sprintf(robot_number, "%d", robot_id_u);
  
   for(i=0; i<FLOCK_SIZE; i++) {
      initialized[i] = 0;       // Set initialization to 0 (= not yet initialized)
   }
  
    printf("Reset: robot %s\n",robot_number);
}


/*
 * Keep given int number within interval {-limit, limit}
 */
void limit(int *number, int limit) {
   if (*number > limit)
      *number = limit;
   if (*number < -limit)
      *number = -limit;
}

/*
 * Updates robot position with wheel speeds
 */
void update_self_motion(int msl, int msr) { 
   float theta = my_position[2];
   //printf("update self motion, %d, %d\n", msl, msr);
   // Compute deltas of the robot
   float dr = (float)msr * SPEED_UNIT_RADS * WHEEL_RADIUS * DELTA_T;
   float dl = (float)msl * SPEED_UNIT_RADS * WHEEL_RADIUS * DELTA_T;
   float du = (dr + dl)/2.0;
   float dtheta = (dr - dl)/AXLE_LENGTH;
  
   // Compute deltas in the environment
   float dx = -du * sinf(theta);
   float dz = -du * cosf(theta);
  
   // Update position
   my_position[0] += dx;
   my_position[1] += dz;
   my_position[2] += dtheta;
  
   // Keep orientation within 0, 2pi
   if (my_position[2] > 2*M_PI) 
      my_position[2] -= 2.0*M_PI;
   if (my_position[2] < 0) 
      my_position[2] += 2.0*M_PI;
   
   //printf("orientation: %f\n", my_position[2]);
}

/*
 * Computes wheel speed given a certain X,Z speed
 */
void compute_wheel_speeds(int *msl, int *msr) 
{
   // Compute wanted position from Reynold's speed and current location
   float x = speed[robot_id][0]*cosf(my_position[2]) - speed[robot_id][1]*sinf(my_position[2]); // x in robot coordinates
   float z = -speed[robot_id][0]*sinf(my_position[2]) - speed[robot_id][1]*cosf(my_position[2]); // z in robot coordinates
   
   float Ku = 0.2;   // Forward control coefficient
   float Kw = 10.0;  // Rotational control coefficient
   float range = sqrtf(x*x + z*z);  // Distance to the wanted position
   float bearing = -atan2(x, z); // Orientation of the wanted position
   
   // Compute forward control
   float u = Ku*range*cosf(bearing);
   // Compute rotational control
   float w = Kw*range*sinf(bearing);
   
   // Convert to wheel speeds!
   *msl = 50*(u - AXLE_LENGTH*w/2.0) / WHEEL_RADIUS;
   *msr = 50*(u + AXLE_LENGTH*w/2.0) / WHEEL_RADIUS;
   limit(msl,MAX_SPEED);
   limit(msr,MAX_SPEED);
}


/*
 *  Update speed according to Reynold's rules
 */

void reynolds_rules() {
   int i, j, k;         // Loop counters
   float rel_avg_loc[2] = {0,0}; // Flock average positions
   float rel_avg_speed[2] = {0,0};  // Flock average speeds
   float cohesion[2] = {0,0};
   float dispersion[2] = {0,0};
   float consistency[2] = {0,0};
   
   /* Compute averages over the whole flock */
   for (j=0;j<2;j++) {
      for(i=0; i<FLOCK_SIZE; i++) 
      {
         // don't consider yourself for the average
         if (i != robot_id) 
         {
            // don't take really old values
            if( (maxTimestamp - timestamp[i]) < MAX_TIMEOUT ){
               rel_avg_speed[j] += relative_speed[i][j];
               rel_avg_loc[j] += relative_pos[i][j];
            }
         }        
      }
   }
   
   
   /* Rule 1 - Aggregation/Cohesion: move towards the center of mass */
    
    for (j=0;j<2;j++) { 
      // If center of mass is too far
        if (fabs(rel_avg_loc[j])> RULE1_THRESHOLD) 
      {     
         cohesion[j] = rel_avg_loc[j] ;  // Relative distance to the center of the swarm
      }  
   }

   /* Rule 2 - Dispersion/Separation: keep far enough from flockmates */
   for (k=0;k<FLOCK_SIZE;k++) {
      if (k != robot_id) {        // Loop on flockmates only
         // If neighbor k is too close
         if (sqrt(pow(relative_pos[k][0],2)+pow(relative_pos[k][1],2)) < RULE2_THRESHOLD) 
         {
            for (j=0;j<2;j++) 
            {
               dispersion[j] -= relative_pos[k][j]; // Relative distance to k
            }
         }
      }
   }
  
   /* Rule 3 - Consistency/Alignment: match the speeds of flockmates */
   for (j=0;j<2;j++) {
      consistency[j] =  rel_avg_speed[j]; // difference speed to the average
    }

    // aggregation of all behaviors with relative influence determined by weights
    for (j=0;j<2;j++) {
      speed[robot_id][j] = cohesion[j] * RULE1_WEIGHT;
        speed[robot_id][j] +=  dispersion[j] * RULE2_WEIGHT;
        speed[robot_id][j] +=  consistency[j] * RULE3_WEIGHT;
        speed[robot_id][j] += (migr[j]-my_position[j]) * MIGRATION_WEIGHT;
     }
}

/*
 * sends a ping of a certain robot value
 * without updating the timestamp
*/
void send_ping_robot(int other_robot_id){
   char out[32];
   //strcpy(out,robot_number);  // in the ping message we send the name of the robot.
   sprintf(out, "%d", other_robot_id);
   sprintf( (out+8), "%u", timestamp[other_robot_id] );
   
   //printf("timestamp sent: %s\n", out+8);
   
   if(other_robot_id != robot_id){
      sprintf( (out+16), "%3.4f", relative_pos[other_robot_id][0] );
      sprintf( (out+24), "%3.4f", relative_pos[other_robot_id][1] ); 
   }
   else{ //no relative pos with the current robot
      sprintf( (out+16), "%3.4f", 0.0f );
      sprintf( (out+24), "%3.4f", 0.0f );
   }
   //printf("%f.4\n", relative_pos[other_robot_id][0] );
   
   wb_emitter_send(emitter2,out,32); 
}

/*
 * send a ping about the content of all the table
*/
void send_ping_all(void){
   //send my position
   timestamp[robot_id] = maxTimestamp+1;
   send_ping_robot(robot_id);
      
   int i = 0;
   for(; i < FLOCK_SIZE; i++){
      if(i != robot_id){
         send_ping_robot(i);
      }
   }
}

/*
 * processing all the received ping messages, and calculate range and bearing to the other robots
 * the range and bearing are measured directly out of message RSSI and direction
*/
void process_received_ping_messages(void)
{
    const double *message_direction;
    double message_rssi; // Received Signal Strength indicator
    double theta;
    double range;
    char *inbuffer;   // Buffer for the receiver node
    int other_robot_id;
    unsigned int recv_timestamp;
   
   while (wb_receiver_get_queue_length(receiver2) > 0) {
      inbuffer = (char*) wb_receiver_get_data(receiver2);
      message_direction = wb_receiver_get_emitter_direction(receiver2);
      message_rssi = wb_receiver_get_signal_strength(receiver2);
      
      //should be x and z position (y is up)
      double x = message_direction[0];
      double z = message_direction[2];
      
      //printf("ROBOT %d: message_direction: %f, %f, %f\n", robot_id, message_direction[0], message_direction[1], message_direction[2]);
      
      theta =  -atan2(z,x);
      theta = theta + my_position[2]; // find the relative theta;
      range = sqrt((1/message_rssi)); 

      //other_robot_id = (int)(inbuffer[5]-'0');  // old: since the name of the sender is in the received message. Note: this does not work for robots having id bigger than 9!
      sscanf(inbuffer, "%d", &other_robot_id);
      sscanf( (inbuffer+8), "%u", &recv_timestamp);
      
      //printf("recieved: id%d, ts%d\n", other_robot_id, recv_timestamp);
      
      //only update the position if there is something new
      if(timestamp[other_robot_id] < recv_timestamp){
         
         if(maxTimestamp < recv_timestamp){
            maxTimestamp = recv_timestamp;
         }
         
         float relativ_pos_x = 0.0f;
         float relativ_pos_z = 0.0f;
         
         sscanf( (inbuffer+16), "%f", &relativ_pos_x);
         sscanf( (inbuffer+24), "%f", &relativ_pos_z);
         
         // Get position update
         prev_relative_pos[other_robot_id][0] = relative_pos[other_robot_id][0];
         prev_relative_pos[other_robot_id][1] = relative_pos[other_robot_id][1];
         
         //get relativ_pos of the sending robot
         relative_pos[other_robot_id][0] = range*cos(theta);  // relative x pos
         relative_pos[other_robot_id][1] = -1.0 * range*sin(theta);   // relative y pos
         
         //add the relativ pos the sending robot gave us
         relative_pos[other_robot_id][0] += relativ_pos_x;
         relative_pos[other_robot_id][1] += relativ_pos_z;
         
         //printf("Robot %s, from robot %d, x: %g, y: %g, theta %g, my theta %g\n",robot_name,other_robot_id,relative_pos[other_robot_id][0],relative_pos[other_robot_id][1],my_position[2]*180.0/3.141592,my_position[2]*180.0/3.141592);
         
         relative_speed[other_robot_id][0] = (1/DELTA_T)*(relative_pos[other_robot_id][0]-prev_relative_pos[other_robot_id][0]);
         relative_speed[other_robot_id][1] = (1/DELTA_T)*(relative_pos[other_robot_id][1]-prev_relative_pos[other_robot_id][1]);
         
         timestamp[other_robot_id] = recv_timestamp;
         
         //sends info to other robots
         send_ping_robot(other_robot_id);
         
      }
      wb_receiver_next_packet(receiver2);
   }
}


// the main function
int main(){ 
   int msl, msr;        // Wheel speeds
   int bmsl, bmsr, sum_sensors;  // Braitenberg parameters
   int i;            // Loop counter
   int distances[NB_SENSORS]; // Array for the distance sensor readings
   int max_sens;        // Store highest sensor value
   
   reset();       // Resetting the robot
   
   for(i=0; i<FLOCK_SIZE; i++) {
      timestamp[i] = 0;
   }

   msl = 0; msr = 0; 
   max_sens = 0; 
   
   // Forever
   for(;;){
      bmsl = 0; bmsr = 0;
      sum_sensors = 0;
      max_sens = 0;
                
      /* Braitenberg */
      for(i=0;i<NB_SENSORS;i++) {
         distances[i]=wb_distance_sensor_get_value(ds[i]); //Read sensor values
         sum_sensors += distances[i]; // Add up sensor values
         max_sens = max_sens>distances[i]?max_sens:distances[i]; // Check if new highest sensor value

         // Weighted sum of distance sensor values for Braitenburg vehicle
         bmsr += braiten_weight[i] * distances[i];
         bmsl += braiten_weight[i+NB_SENSORS] * distances[i];
        }

      // Adapt Braitenberg values (empirical tests)
      bmsl/=MIN_SENS; bmsr/=MIN_SENS;
      bmsl+=66; bmsr+=72;
      
      /* Send and get information */
      timestamp[robot_id]++;
      send_ping_robot(robot_id); 
      
      process_received_ping_messages();
      
      //printf("ROBOT %d: wheels %d, %d\n", robot_id, bmsl, bmsr);
               
      // Compute self position
      prev_my_position[0] = my_position[0];
      prev_my_position[1] = my_position[1];
      
      update_self_motion(msl,msr);
      
      speed[robot_id][0] = (1/DELTA_T)*(my_position[0]-prev_my_position[0]);
      speed[robot_id][1] = (1/DELTA_T)*(my_position[1]-prev_my_position[1]);
      
      // Reynold's rules with all previous info (updates the speed[][] table)
      reynolds_rules();
      
      //printf("ROBOT %d: wanted position: %f, %f\n", robot_id, speed[robot_id][0], speed[robot_id][1]);
      
      // Compute wheels speed from reynold's speed
      compute_wheel_speeds(&msl, &msr);
      
      //printf("wheels: %d, %d\n", msl, msr);
      
      // Adapt speed instinct to distance sensor values
      if (sum_sensors > NB_SENSORS*MIN_SENS) {
         msl -= msl*max_sens/(2*MAX_SENS);
         msr -= msr*max_sens/(2*MAX_SENS);
      }
    
      // Add Braitenberg
      msl += bmsl;
      msr += bmsr;
                  
      // Set speed
      wb_differential_wheels_set_speed(msl,msr);
    
      // Continue one step
      wb_robot_step(TIME_STEP);
   }
}  
