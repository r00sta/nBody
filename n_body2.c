/* #### n_body2 #### */
/* An N-body  simulation program using the leapfrog integration method.
 * Data is output to STOUT in the form of the data array which is form-
 * ed as data[timestep][particle #][(t,x,y,z,V,K)]. * 
 * */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>


/* #### Define Program Parameters #### */

#define NumP 10 														/*Number of Particles*/
#define dt 0.1 															/*Timestep in days*/
#define Ndt 100000 														/*Number of Timesteps*/
#define G 6.67E-11 														/*Gravitational Constant*/
#define e 1E11 															/*Epsilon Value*/
#define prec 100 														/*Set Output Precision: 1-Full Precision, >1-Less Precise*/
#define size 1E12 														/*Universe Size*/


/* #### Inititalize working arrays #### */

double masses[NumP]={0.0}; 												/*Particle Masses in Solar Mass Units*/
double pos[NumP][3]={{0.0}}; 											/*Particle Positions in m (x,y,z)*/
double vel[NumP][3]={{0.0}}; 											/*Particle Velocities in ms^-1 (v_x,v_y,v_z)*/
double new_pos[NumP][3]={{0.0}}; 										/*New Particle Positions in m (x,y,z)*/
double new_vel[NumP][3]={{0.0}}; 										/*New Particle Velocities in ms^-1 (v_x,v_y,v_z)*/
double time[NumP]={0.0}; 												/*Time Array in Days*/
double data[Ndt][NumP][6]={{{0.0}}}; 									/*Output Data array [timestep][particle number][x,y,z,t,Ek,Ep]*/
double dpos[NumP][4]={{0.0}}; 											/*Distance between particles for a given particle [particle #][dx,dy,dz,dr]*/
double force[NumP][3]={{0.0}}; 											/*Force vectors [i][Fx,Fy,Fz]*/
double Fx, Fy, Fz;														/*Force vectors*/
double kin_en[NumP]={0.0}; 												/*Particles Kinetic Energy*/
double pot_en[NumP]={0.0}; 												/*Particles Potential Energy*/


/* #### Functions #### */

int Iterate(particle_i,timestep)										/*Particle iteration function*/
{
	Fx = Fy = Fz = 0.0;													/*Initiate forces to zero*/
	
	pot_en[particle_i] = 0.0;											/*Set potential energy for particle_i to zero*/
	
	int i;
	
	for (i = 0; i < NumP; i++)											/*Loop over particles*/
	{
		
		int j;
		
		for (j=0; j < 3; j++)											/*Loop over x,y,z*/
		{
			dpos[i][j] = pos[particle_i][j] - pos[i][j];				/*Caclulate particle seperation along axis*/
		}
		
		/*Calculate particle radial seperation*/
		dpos[i][3] = sqrt( pow(dpos[i][0],2) + pow(dpos[i][1],2) + pow(dpos[i][2],2) );
		
		/*Calculate Force in x-direction and append to Fx*/
		force[particle_i][0] = - (G * masses[particle_i] * masses[i] * dpos[i][0]) / pow( sqrt(pow(dpos[i][3],2) + pow(e,2)) ,3);
		Fx = Fx + force[particle_i][0];
		
		/*Calculate Force in y-direction and append to Fy*/
		force[particle_i][1] = - (G * masses[particle_i] * masses[i] * dpos[i][1]) / pow( sqrt(pow(dpos[i][3],2) + pow(e,2)) ,3);
		Fy = Fy + force[particle_i][1];
		
		/*Calculate Force in z-direction and append to Fz*/
		force[particle_i][2] = - (G * masses[particle_i] * masses[i] * dpos[i][2]) / pow( sqrt(pow(dpos[i][3],2) + pow(e,2)) ,3);
		Fz = Fz + force[particle_i][2];				
		
		if (dpos[i][3] != 0.0) 											/*Skip self potetials*/
		{
			/*Calculate potential energy of particle_i*/
			pot_en[particle_i] = pot_en[particle_i] - (G * masses[particle_i] * masses[i])/(dpos[i][3]+1);
		}		
	}
	
	data[timestep][particle_i][5] = pot_en[particle_i];					/*Add potetial energy to data array*/
	
	/*Calculate new velocity in each direction and add to new_vel array*/
	new_vel[particle_i][0] = vel[particle_i][0] + ((Fx * dt * 3600 * 24)/(masses[particle_i]));
	new_vel[particle_i][1] = vel[particle_i][1] + ((Fy * dt * 3600 * 24)/(masses[particle_i]));
	new_vel[particle_i][2] = vel[particle_i][2] + ((Fz * dt * 3600 * 24)/(masses[particle_i]));
	
	/*Calculate new x positon and add to data array*/
	new_pos[particle_i][0] = pos[particle_i][0] + (new_vel[particle_i][0] * dt * 3600 * 24);
	data[timestep][particle_i][0] = new_pos[particle_i][0];
	
	/*Calculate new y positon and add to data array*/
	new_pos[particle_i][1] = pos[particle_i][1] + (new_vel[particle_i][1] * dt * 3600 * 24);
	data[timestep][particle_i][1] = new_pos[particle_i][1];
	
	/*Calculate new z positon and add to data array*/
	new_pos[particle_i][2] = pos[particle_i][2] + (new_vel[particle_i][2] * dt * 3600 * 24);
	data[timestep][particle_i][2] = new_pos[particle_i][2];
	
	/*Calculate velocities half a timestep for kinetic energy calcs*/
	double vel_x = vel[particle_i][0] + ((Fx * dt * 3600 * 24 * 0.5)/(masses[particle_i]));
	double vel_y = vel[particle_i][1] + ((Fy * dt * 3600 * 24 * 0.5)/(masses[particle_i]));
	double vel_z = vel[particle_i][2] + ((Fz * dt * 3600 * 24 * 0.5)/(masses[particle_i]));
	double vel_abs = sqrt(pow(vel_x,2) + pow(vel_y,2) + pow(vel_z,2));
	
	/*Calculate kinectic energy and add to kin_en array*/
	kin_en[particle_i] = masses[particle_i]*pow(vel_abs,2);
	data[timestep][particle_i][4] = kin_en[particle_i];
	
	return 0;
}	

double rand_range(double min_n, double max_n)
{
    return (double)rand()/RAND_MAX * (max_n - min_n) + min_n;
}

int RandomPos()
{
	double max = 1.0;
	double min = -1.0;
	int particles = 0;
	for (particles=0; particles < NumP; particles++)
	{
		masses[particles] = 2E27;
		pos[particles][0] = rand_range(min,max)*1.0E12;
		pos[particles][1] = rand_range(min,max)*1.0E12;
		pos[particles][2] = rand_range(min,max)*1.0E12;
	}
	return 0;
}


/* #### Main Program #### */

int main()
{
	RandomPos(); 														/*Asign random positions to particles*/
	
   	int t; 																/*Initiate timestep counter*/
   	for (t = 0; t < Ndt; t++) 											/*Loop through timesteps*/
   	{
		
    	int particle_i; 												/*Initiate particle counter*/
    	
    	for (particle_i = 0; particle_i < NumP; particle_i++) 			/*Loop through particles*/
    	{
			data[t][particle_i][3] = t * dt; 							/*Add time to data array*/
			Iterate(particle_i,t); 										/*Iterate particle_i*/
		}
		
		int i;
		
		for (i=0; i < NumP; i++) 										/*Iterate over particles*/
		{
			int j;
			for (j=0; j < 3; j++)
			{
				pos[i][j] = new_pos[i][j]; 								/*Asign new positions to position array*/
				vel[i][j] = new_vel[i][j]; 								/*Asign new velocities to velocity array*/
			}
		}
   	}
   	
   	/*Data Output to stdio*/
   	
   	printf ("NumP,dt,Ndt,Precision, ,\n");								
   	printf ("%i,%f,%i,%i,%e,\n", NumP, dt, Ndt, prec, size);			/*System Parameters used in analysis*/
   	printf("x,y,z,t,Ek,Ep\n");
   	int i;
   	
   	for (i=0; i < Ndt; i+=prec)											/*Loop through timesteps with specified precision*/
   	{
		
		int j;
		
		for (j=0; j < NumP; j++)										/*Loop over particles*/
		{
			/*Print particle data (x,y,z,t,Ek,Ep)*/
			printf("%e,%e,%e,%e,%e,%e\n" , data[i][j][0], data[i][j][1], data[i][j][2], data[i][j][3], data[i][j][4], data[i][j][5]);
		}   	
   	}
   	   	   	
   	return 0;
}



