#include <stdio.h>
#include <math.h>
#include <string.h>

#include <webots/robot.h>

#include <webots/differential_wheels.h>
#include <webots/distance_sensor.h>
#include <webots/emitter.h>
#include <webots/receiver.h>
#include <webots/compass.h>

//Braitenberg parameters for obstacle avoidance
#define NB_SENSORS		8	  // Number of distance sensors
#define MIN_SENS        350       // Minimum sensibility value
#define MAX_SENS        4096      // Maximum sensibility value
#define MAX_SPEED       800       // Maximum speed
#define FLOCK_SIZE		5	  // Size of flock
#define TIME_STEP		64	  // [ms] Length of time step

#define AXLE_LENGTH 	0.052	  // Distance between wheels of robot (meters)
#define SPEED_UNIT_RADS	0.00628	  // Conversion factor from speed unit to radian per second
#define WHEEL_RADIUS	0.0205	  // Wheel radius (meters)
#define DELTA_T			0.064	  // Timestep (seconds)

#define KP_COEFFICIENT 0.5f
#define ODES_COEFFICIENT 400
#define C_COEFFICIENT 4096.0
#define BETA_COEFFICIENT 3
 
int e_puck_matrix[16] = {17,29,34,10,8,-38,-56,-76,-72,-58,-36,8,10,36,28,18}; // for obstacle avoidance
float sensor_degrees[] = {10.0, 20.0, 45.0, 165.0, 195.0, 270.0, 340.0, 350.0};

WbDeviceTag ds[NB_SENSORS];	// Handle for the infrared distance sensors
WbDeviceTag receiver2;		// Handle for the receiver node
WbDeviceTag emitter2;		// Handle for the emitter node
WbDeviceTag compass2;

int robot_id_u, robot_id;	// Unique and normalized (between 0 and FLOCK_SIZE-1) robot ID

float relative_pos[FLOCK_SIZE][3];	// relative X, Z, Theta of all robots
float neighboors_bearing[FLOCK_SIZE]; // Bearing of the neighboors
float prev_relative_pos[FLOCK_SIZE][3];	// Previous relative  X, Z, Theta values
float my_position[3];     		// X, Z, Theta of the current robot
float prev_my_position[3];  		// X, Z, Theta of the current robot in the previous time step
float speed[FLOCK_SIZE][2];		// Speeds calculated with Reynold's rules
float relative_speed[FLOCK_SIZE][2];	// Speeds calculated with Reynold's rules
int initialized[FLOCK_SIZE];		// != 0 if initial positions have been received
float migr[2] = {0,-50};	        // Migration vector
char* robot_name; 


/*
 * Reset the robot's devices and get its ID
 */
static void reset() 
{
	wb_robot_init();

	receiver2 = wb_robot_get_device("receiver2");
	emitter2 = wb_robot_get_device("emitter2");
	compass2 = wb_robot_get_device("compass2");
	
	int i;
	char s[4]="ps0";
	for(i=0; i<NB_SENSORS;i++) {
		ds[i]=wb_robot_get_device(s);	// the device name is specified in the world file
		s[2]++;				// increases the device number
	}
	robot_name=(char*) wb_robot_get_name(); 

	for(i=0;i<NB_SENSORS;i++)
    	wb_distance_sensor_enable(ds[i],64);

	wb_receiver_enable(receiver2,64);
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
	double rad = atan2(migr[0], migr[1]);
	double bearing = (rad - 1.5708) / M_PI * 180.0;
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
	float h[2];
	float a[2], ac[2];
	float p[2];

	float u = 0;
	float norm = 0;
	float dot = 0;
	float w = 0;
	float fk = 0;

	// Heading vector H
	h[0] = h[1] = 0.0f;
	for(i = 0; i < FLOCK_SIZE; ++i) {
		if (i != robot_id) {
			h[0] += cosf((get_bearing_in_degrees() - neighboors_bearing[i]) * M_PI / 180.0);
			h[1] += sinf((get_bearing_in_degrees() - neighboors_bearing[i]) * M_PI / 180.0);
		}
	}
	norm = sqrtf(h[0]*h[0] + h[1]*h[1]);
	h[0] /= norm;
	h[1] /= norm;
	printf("H : %.4lf %.4lf\n", h[0], h[1]);	

	// Proximal control behavior fk and P
	p[0] = p[1] = .0f;
	for(i = 0; i < NB_SENSORS; ++i) {
		fk = ds[i]-ODES_COEFFICIENT;
		fk = (fk*fk)/C_COEFFICIENT;

		if (ds[i] >= ODES_COEFFICIENT)
			fk *= -1.0f;

		p[0] += cosf(sensor_degrees[i] * M_PI / 180.0)*fk;
		p[1] += sinf(sensor_degrees[i] * M_PI / 180.0)*fk;
	}
	p[0] /= ((float)NB_SENSORS);
	p[1] /= ((float)NB_SENSORS);
	printf("P : %.4lf %.4lf\n", p[0], p[1]);

	// Compute final desired vector a
	a[0] = a[1] = 0.0f;
	a[0] = h[0] + BETA_COEFFICIENT * p[0];
	a[1] = h[1] + BETA_COEFFICIENT * p[1];
	norm = sqrtf(a[0]*a[0] + a[1]*a[1]); 
	a[0] /= norm;
	a[1] /= norm;
	//printf("A : %.4lf %.4lf\n", a[0], a[1]);

	// Motion control U

	// Current heading
	ac[0] = speed[robot_id][0];
	ac[1] = speed[robot_id][1];
	norm = sqrtf(ac[0]*ac[0] + ac[1]*ac[1]);
	ac[0] /= norm;
	ac[1] /= norm;

	dot = ac[0]*a[0] + ac[1]*a[1];
	u = (dot >= 0.0) ? dot * 0.1287 : 0.0;
	//printf("U : %.4lf\n", u);
	// Angular velocity w

	w = KP_COEFFICIENT * (atanf(ac[1] / ac[0]) - atanf(a[1] / a[0]));
	//printf("W : %.4lf %.4lf %.4lf\n", w, a[0], a[1]);

	// Compute wheel speed
	*msr = ((int)(50*(u - AXLE_LENGTH*w/2.0) / WHEEL_RADIUS));
	*msl = ((int)(50*(u + AXLE_LENGTH*w/2.0) / WHEEL_RADIUS));

	//*msr *= 1000.0 / (2.0 * M_PI);
	//*msl *= 1000.0 / (2.0 * M_PI);
	limit(msl,MAX_SPEED);
	limit(msr,MAX_SPEED);

	printf("MSL/R %d %d\n", *msl, *msr);
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
	wb_emitter_send(emitter2,out,strlen(out)+1); 
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
	
	while (wb_receiver_get_queue_length(receiver2) > 0) {
		inbuffer = (char*) wb_receiver_get_data(receiver2);
		message_direction = wb_receiver_get_emitter_direction(receiver2);
		message_rssi = wb_receiver_get_signal_strength(receiver2);
		double y = message_direction[2];
		double x = message_direction[1];

        theta =	-atan2(y,x);
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

		relative_speed[other_robot_id][0] = (1/DELTA_T)*(relative_pos[other_robot_id][0]-prev_relative_pos[other_robot_id][0]);
		relative_speed[other_robot_id][1] = (1/DELTA_T)*(relative_pos[other_robot_id][1]-prev_relative_pos[other_robot_id][1]);		
		wb_receiver_next_packet(receiver2);
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
/*                
		// Braitenberg
		for(i=0;i<NB_SENSORS;i++) {
			distances[i]=wb_distance_sensor_get_value(ds[i]); //Read sensor values
			sum_sensors += distances[i]; // Add up sensor values
			max_sens = max_sens>distances[i]?max_sens:distances[i]; // Check if new highest sensor value

			// Weighted sum of distance sensor values for Braitenburg vehicle
			bmsr += e_puck_matrix[i] * distances[i];
			bmsl += e_puck_matrix[i+NB_SENSORS] * distances[i];
        }

		// Adapt Braitenberg values (empirical tests)
        bmsl/=MIN_SENS; bmsr/=MIN_SENS;
		// bmsl+=66; bmsr+=72;
*/              
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
    
		// Adapt speed instinct to distance sensor values
		if (sum_sensors > NB_SENSORS*MIN_SENS) {
			msl -= msl*max_sens/(2*MAX_SENS);
			msr -= msr*max_sens/(2*MAX_SENS);
		}
    
		// Add Braitenberg
		//msl += bmsl;
		//msr += bmsr;
                  
		// Set speed
		wb_differential_wheels_set_speed(msl,msr);
    
		// Continue one step
		wb_robot_step(TIME_STEP);
	}
}  
  
