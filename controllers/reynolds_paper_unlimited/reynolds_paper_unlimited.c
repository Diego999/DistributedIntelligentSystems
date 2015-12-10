#include <stdio.h>
#include <math.h>
#include <string.h>

#include <webots/robot.h>

#include <webots/differential_wheels.h>
#include <webots/distance_sensor.h>
#include <webots/emitter.h>
#include <webots/receiver.h>
#include <webots/compass.h>

/*
 * SIMULATION PARAMETERS
 */

#define FLOCK_SIZE		5	     // Size of flock
#define TIME_STEP		64	     // [ms] Length of time step

/*
 * E-PUCK PARAMETERS
 */
#define MIN_SENS        350      // Minimum sensibility value
#define MAX_SENS        4096     // Maximum sensibility value 
#define NB_SENSORS		8	     // Number of distance sensors
#define AXLE_LENGTH 	0.052	 // Distance between wheels of robot (meters)
#define SPEED_UNIT_RADS	0.00628	 // Conversion factor from speed unit to radian per second
#define WHEEL_RADIUS	0.0205	 // Wheel radius (meters)
#define DELTA_T			0.064    // Timestep (seconds)
#define MAX_SPEED_MS    0.1287   // Maximum linear speed for the e-pucks
#define MAX_SPEED       800      // Maximum speed in tics
float sensor_degrees[] = {-18.435, -48.652, -90.000, -116.565, +116.565, +90.000, +48.652, +18.435}; // position in degrees for sensor position

/*
 * CONTROL PARAMETERS
 */
#define COEF_KP 1.0    // proportional gain for the angular velocity controller
#define COEF_B  8.0    // gain on the virtual force
#define COEF_R  1.0    // gain on the migration urge
#define COEF_S  5.0    // gain on the cohesion
#define COEF_T  0.1    // Threshold
#define DISPLAY 1       // Display on/off for e-puck state values while controled


/*
 * MIGRATION URGE 
 */
float migr[2] = {+100.0,0.0}; // Migration vector

/*
 * OTHER PARAMETERS
 */
WbDeviceTag ds[NB_SENSORS];	 // Handle for the infrared distance sensors
WbDeviceTag receiver;       // Handle for the receiver nCOEF_oDe
WbDeviceTag emitter;		 // Handle for the emitter nCOEF_oDe
WbDeviceTag compass2;        // Handle for the compass

int robot_id_u, robot_id;	// Unique and normalized (between 0 and FLOCK_SIZE-1) robot ID

float relative_pos[FLOCK_SIZE][3];	    // relative X, Z, Theta of all robots
float neighboors_bearing[FLOCK_SIZE];   // Bearing of the neighboors
float prev_relative_pos[FLOCK_SIZE][3];	// Previous relative  X, Z, Theta values
float my_position[3];     		        // X, Z, Theta of the current robot
float prev_my_position[3];  		    // X, Z, Theta of the current robot in the previous time step
float speed[FLOCK_SIZE][2];		        // Speeds calculated with Reynold's rules
float relative_speed[FLOCK_SIZE][2];	// Speeds calculated with Reynold's rules
int initialized[FLOCK_SIZE];		    // != 0 if initial positions have been received
char* robot_name; 


/*
 * Reset the robot's devices and get its ID
 */
static void reset() 
{
	wb_robot_init();

	receiver = wb_robot_get_device("receiver");
	emitter  = wb_robot_get_device("emitter");
	compass2  = wb_robot_get_device("compass2");
	
	int i;
	char s[4]="ps0";
	for(i=0; i<NB_SENSORS;i++) {
		ds[i]=wb_robot_get_device(s);  // the device name is specified in the world file
		s[2]++;				           // increases the device number
	}
	robot_name=(char*) wb_robot_get_name(); 

	for(i=0;i<NB_SENSORS;i++)
    	wb_distance_sensor_enable(ds[i],64);

	wb_receiver_enable(receiver,64);
	wb_compass_enable(compass2, 64);
	
	//Reading the robot's name. Pay attention to name specification when adding robots to the simulation!
	sscanf(robot_name,"epuck%d",&robot_id_u); // read robot id from the robot's name
	robot_id = robot_id_u%FLOCK_SIZE;	  // normalize between 0 and FLOCK_SIZE-1
  
	for(i=0; i<FLOCK_SIZE; i++) {
		initialized[i] = 0;		  // Set initialization to 0 (= not yet initialized)
	}
  
    printf("Reset: robot %d\n",robot_id_u);
}

double get_bearing_in_degrees() {
	const double *north = wb_compass_get_values(compass2);
	double rad = atan2(north[0], north[2]);
	double bearing = (rad - 1.5708) / M_PI * 180.0;
	if (bearing < 0.0)
		bearing = bearing + 360.0;
	return bearing;
}

double get_migr_bearing_in_degrees() {
	//const double *north = wb_compass_get_values(compass2);
	double rad = atan2(migr[1], migr[0]);
	double bearing = (rad) / M_PI * 180.0;
	if (bearing < 0.0)
		bearing = bearing + 360.0;
	return bearing;
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
}

void flocking_behavior(int *msl, int *msr) 
{
	int i;
	
	float a[2]; // desired heading vector 
	float b[2]; // current heading vector
	float h[2]; // haeding alignment vector
	float p[2]; // proximal control vector
	float r[2]; // migration control vector
	float s[2]; // flocking control vector

	float u  = 0.0; // linear velocity of the e-puck
	float w  = 0.0; // angular velocity of the e-puck
	float Nl = 0.0; // angular velocity of the left wheel
	float Nr = 0.0; // angular velocity of the right wheel
	
	float fk   = 0.0; // buffer variable to compute the virtual force from IR sensor
	float norm = 0.0; // buffer variable for vector normalization
	float sens = 0.0; // buffer variable to get senser values
	float xloc = 0.0; // buffer variable for the robots relative x position
	float zloc = 0.0; // buffer variable for the robots relative z position
	float aloc = 0.0; // buffer variable for the robots relative anlge
	float dloc = 0.0; // buffer variable for the robots relative distance

	/* 
	 * HEADING VECTOR
	 */ 
	h[0] = h[1] = 0.0f;
	for(i = 0; i < FLOCK_SIZE; ++i) {
		if (i != robot_id) {
			h[0] += cosf( ( get_bearing_in_degrees() - neighboors_bearing[i] ) * M_PI / 180.0 );
			h[1] += sinf( ( get_bearing_in_degrees() - neighboors_bearing[i] ) * M_PI / 180.0 );
		}
	}

	norm = sqrtf(h[0]*h[0] + h[1]*h[1]);
	h[0] /= norm;
	h[1] /= norm;

	/*
	 * PROXIMAL CONTROL
	 */
	p[0] = p[1] = .0f;
	for(i = 0; i < NB_SENSORS; ++i) {
		sens = wb_distance_sensor_get_value(ds[i]);
		if (sens >= MIN_SENS) fk = -(sens*sens)/MAX_SENS;
    	if (sens <  MIN_SENS) fk = 0.0;
		p[0] += cosf( sensor_degrees[i] * M_PI / 180.0 ) * fk / ((float)NB_SENSORS);
		p[1] += sinf( sensor_degrees[i] * M_PI / 180.0 ) * fk / ((float)NB_SENSORS);
	}

	s[0] = s[1] = .0f;
	for(i = 0; i < FLOCK_SIZE; i++) {
		fk = 0.0;
		xloc = -relative_pos[i][1];
		zloc = -relative_pos[i][0];
		dloc = sqrtf(xloc*xloc + zloc*zloc);
		if (i != robot_id)
		{
			if ( xloc >= 0.0) aloc = +atanf(zloc/xloc) * 180/M_PI;
			if ( xloc <  0.0) aloc = +atanf(zloc/xloc) * 180/M_PI + 180;
			if ( dloc >= COEF_T)  fk = +(dloc - COEF_T)*(dloc - COEF_T);
			if ( dloc <  COEF_T)  fk = -(dloc - COEF_T)*(dloc - COEF_T); 
			s[0] += cosf(aloc * M_PI/180) * fk;
			s[1] += sinf(aloc * M_PI/180) * fk;
		}		
	}
	
	/*
	 * MIGRATION URGE
	 */
  	r[0] = cosf(( get_bearing_in_degrees() - get_migr_bearing_in_degrees() - 90)* M_PI/180);
  	r[1] = sinf(( get_bearing_in_degrees() - get_migr_bearing_in_degrees() - 90)* M_PI/180);
	/*
	 * HEADING VECTOR (desired = a, current = b)
	 */
	a[0] = h[0] + COEF_B * p[0] + COEF_R * r[0] + COEF_S* s[0];
	a[1] = h[1] + COEF_B * p[1] + COEF_R * r[1] + COEF_S* s[1];
	norm = sqrtf(a[0]*a[0] + a[1]*a[1]);
	a[0] /= norm;
	a[1] /= norm;
	b[0] = cosf( 0.0 * M_PI / 180 );
	b[1] = sinf( 0.0 * M_PI / 180 );

	/*
	 * MOTION CONTROL
	 */
	float phase_a = atanf(a[1] / a[0]);
	float phase_b = atanf(b[1] / b[0]);
	float dphase  = phase_b - phase_a;

	u = (a[0]*b[0] + a[1]*b[1]) * MAX_SPEED_MS;
	w = dphase * COEF_KP ;
	if (u < 0.0){
		u = 0.0;
		if (phase_a >  0) dphase = phase_b - (phase_a - M_PI);
		if (phase_a <= 0) dphase = phase_b - (phase_a + M_PI);
		w = dphase * COEF_KP;
	}

	/*
	 * WHEEL VELOCITY
	 */
	Nl   = (u + AXLE_LENGTH * w / 2.0 ) / (2 * M_PI * WHEEL_RADIUS);
	Nr   = (u - AXLE_LENGTH * w / 2.0 ) / (2 * M_PI * WHEEL_RADIUS);
	*msl = (int) ( Nl / ( SPEED_UNIT_RADS )) ; 
	*msr = (int) ( Nr / ( SPEED_UNIT_RADS )) ;
	limit(msl,MAX_SPEED);
	limit(msr,MAX_SPEED);

	/*
	 * PRINT
	 */
	if (DISPLAY == 2)
	{
		if (robot_id == 3){
		printf(" ------ ROBOT ID ------ \n");
		printf(" robot id : %d \n " , robot_id);
  		printf(" ------ VECTORS ------ \n ");
		printf(" h        : %.6lf , %.6lf \n" , h[0], h[1]);	
		printf(" p        : %.6lf , %.6lf \n" , p[0], p[1]);
		printf(" r        : %.6lf , %.6lf \n" , r[0], r[1]);	
		printf(" a        : %.6lf , %.6lf \n" , a[0], a[1]);
		printf(" b        : %.6lf , %.6lf \n" , b[0], b[1]);
		printf(" s        : %.6lf , %.6lf \n" , s[0], s[1]);
		printf(" ------ CONTROL ------ \n ");
		printf(" pa , pb  : %.6lf , %.6lf \n" , phase_a*(180/M_PI), phase_b*(180/M_PI));
		printf(" u  , w   : %.6lf , %.6lf \n" , u, w);
		printf(" Nl , Nr  : %.6lf , %.6lf \n" , Nl , Nr);
 		printf(" MSL/R    : %.6d  , %.6d  \n", *msl, *msr);
		printf(" --------------------- \n "); 
		}
	}
}

/*
 *  each robot sends a ping message, so the other robots can measure relative range and bearing to the sender.
 *  the message contains the robot's name
 *  the range and bearing will be measured directly out of message RSSI and direction
*/
void send_ping(void)  
{
    char out[100];
	sprintf(out,"%s;%.5f", robot_name, get_bearing_in_degrees());
	wb_emitter_send(emitter,out,strlen(out)+1); 
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
	char *inbuffer;	// Buffer for the receiver node
    int other_robot_id;
    double other_robot_bearing;
	
	while (wb_receiver_get_queue_length(receiver) > 0) {
		inbuffer = (char*) wb_receiver_get_data(receiver);
		message_direction = wb_receiver_get_emitter_direction(receiver);
		message_rssi = wb_receiver_get_signal_strength(receiver);
		
		//should be x and z position (y is up)
		double x = message_direction[0];
		double z = message_direction[2];
		
      
      //printf("ROBOT %d: message_direction: %f, %f, %f\n", robot_id, message_direction[0], message_direction[1], message_direction[2]);
      
        theta =	-atan2(z,x);
        theta = theta + my_position[2]; // find the relative theta;
		range = sqrt((1/message_rssi)); 

		other_robot_id = (int)(inbuffer[5]-'0');  // since the name of the sender is in the received message. Note: this does not work for robots having id bigger than 9!
		other_robot_bearing = 0;
		sscanf(inbuffer+7, "%lf", &other_robot_bearing);

		// Get position update
		prev_relative_pos[other_robot_id][0] = relative_pos[other_robot_id][0];
		prev_relative_pos[other_robot_id][1] = relative_pos[other_robot_id][1];

		relative_pos[other_robot_id][0] = range*cos(theta);  // relative x pos
		relative_pos[other_robot_id][1] = -1.0 * range*sin(theta);   // relative y pos

		// Get bearing
		neighboors_bearing[other_robot_id] = other_robot_bearing;

		//printf("Robot %s, from robot %d, x: %g, y: %g, theta %g, my theta %g\n",robot_name,other_robot_id,relative_pos[other_robot_id][0],relative_pos[other_robot_id][1],my_position[2]*180.0/3.141592,my_position[2]*180.0/3.141592);

		relative_speed[other_robot_id][0] = (1/DELTA_T)*(relative_pos[other_robot_id][0]-prev_relative_pos[other_robot_id][0]);
		relative_speed[other_robot_id][1] = (1/DELTA_T)*(relative_pos[other_robot_id][1]-prev_relative_pos[other_robot_id][1]);		
		wb_receiver_next_packet(receiver);
	}
}

// the main function
int main(){ 
	int msl, msr;			// Wheel speeds
	int bmsl, bmsr, sum_sensors;	// Braitenberg parameters
	//int i;				// Loop counter
	//int distances[NB_SENSORS];	// Array for the distance sensor readings
	int max_sens;			// Store highest sensor value
	
 	reset();			// Resetting the robot
  
	msl = 0; msr = 0; 
	max_sens = 0; 
	
	// Forever
	for(;;){
        bmsl = 0; bmsr = 0;
        sum_sensors = 0;
		max_sens = 0;
        
		/* Send and get information */
		send_ping();  // sending a ping to other robot, so they can measure their distance to this robot
		
		process_received_ping_messages();
					
		// Compute self position
		prev_my_position[0] = my_position[0];
		prev_my_position[1] = my_position[1];
		
		update_self_motion(msl,msr);
		
		speed[robot_id][0] = (1/DELTA_T)*(my_position[0]-prev_my_position[0]);
		speed[robot_id][1] = (1/DELTA_T)*(my_position[1]-prev_my_position[1]);
    
		// Flocking behavior of the paper with the wheels speed
		flocking_behavior(&msl, &msr);
    
		// Set speed
		wb_differential_wheels_set_speed(msl,msr);
    
		// Continue one step
		wb_robot_step(TIME_STEP);
	}
}  
  
