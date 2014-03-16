/* #### n_body2 #### */
/* An N-body  simulation program using the leapfrog integration method.
 * Data is output to STOUT in the form of the data array which is form-
 * ed as data[timestep][particle #][(t,x,y,z,V,K)]. * 
 * */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h> 


/* #### Define Program Parameters #### */

#define NumP 2 														/*Number of Particles*/
#define dt 1.0													/*Timestep in days*/
#define Ndt 200000														/*Number of Timesteps*/
#define G 6.67E-11 														/*Gravitational Constant*/
#define e 1E10 															/*Epsilon Value*/
#define prec 2000												/*Set Output Precision: 1-Full Precision, >1-Less Precise*/
#define size 1.0E13 												/*Universe Size*/
#define a 0.0															/*a(t) expansion factor. Start value a(t=0).*/
#define M 2.0E28														/*Particle Mass*/


/* #### Inititalize working arrays #### */

double masses[NumP]={0.0}; 												/*Particle Masses in Solar Mass Units*/
double pos[NumP][3]={{0.0}}; 											/*Particle Positions in m (x,y,z)*/
double vel[NumP][3]={{0.0}}; 											/*Particle Velocities in ms^-1 (v_x,v_y,v_z)*/
double new_pos[NumP][3]={{0.0}}; 										/*New Particle Positions in m (x,y,z)*/
double new_vel[NumP][3]={{0.0}}; 										/*New Particle Velocities in ms^-1 (v_x,v_y,v_z)*/
double time[NumP]={0.0}; 												/*Time Array in Days*/
double data[Ndt][NumP][6]={{{0.0}}}; 									/*Output Data array [timestep][particle number][x,y,z,t,Ek,Ep]*/
double dpos[NumP][NumP][4]={{{0.0}}}; 											/*Distance between particles for a given particle [particle #][dx,dy,dz,dr]*/
double force[NumP][NumP][3]={{{0.0}}}; 											/*Force vectors [i][Fx,Fy,Fz]*/
double Fx[NumP] = {0.0};
double Fy[NumP] = {0.0};
double Fz[NumP] = {0.0};														/*Force vectors*/
double kin_en[NumP]={0.0}; 												/*Particles Kinetic Energy*/
double pot_en[NumP]={0.0}; 												/*Particles Potential Energy*/


/* #### Functions #### */


int Calc()
{
	
	int z;
	for (z=0; z < NumP; z++)
	{
		Fx[z] = 0.0;
		Fy[z] = 0.0;
		Fz[z] = 0.0;
		pot_en[z] = 0.0;
	}	
	int i;
	for (i=0; i<(NumP-1); i++)
	{
		int j;
		/*pot_en[i] = 0.0;*/
		for (j=(i+1); j<NumP; j++)
		{
			/*pot_en[j] = 0.0;*/
			
			int k;
			for (k=0; k<3; k++)
			{
				dpos[i][j][k] = pos[i][k] - pos[j][k];
				dpos[j][i][k] = - dpos[i][j][k];
			}
			dpos[i][j][3] = sqrt( pow(dpos[i][j][0],2) + pow(dpos[i][j][1],2) + pow(dpos[i][j][2],2) );
			dpos[j][i][3] = dpos[i][j][3];
			
			/*Calculate Force in x-direction and append to Fx*/
			force[i][j][0] = - (G * masses[i] * masses[j] * dpos[i][j][0]) / pow( sqrt(pow(dpos[i][j][3],2) + pow(e,2)) ,3);
			force[j][i][0] = - force[i][j][0];
			Fx[i] = Fx[i] + force[i][j][0];
			Fx[j] = Fx[j] + force[j][i][0];
			/*Calculate Force in y-direction and append to Fy*/
			force[i][j][1] = - (G * masses[i] * masses[j] * dpos[i][j][1]) / pow( sqrt(pow(dpos[i][j][3],2) + pow(e,2)) ,3);
			force[j][i][1] = - force[i][j][1];
			Fy[i] = Fy[i] + force[i][j][1];
			Fy[j] = Fy[j] + force[j][i][1];
			/*Calculate Force in z-direction and append to Fz*/
			force[i][j][2] = - (G * masses[i] * masses[j] * dpos[i][j][2]) / pow( sqrt(pow(dpos[i][j][3],2) + pow(e,2)) ,3);
			force[j][i][2] = - force[i][j][2];
			Fz[i] = Fz[i] + force[i][j][2];
			Fz[j] = Fz[j] + force[j][i][2];
			pot_en[i] = pot_en[i] - (G * masses[i] * masses[j])/(dpos[i][j][3]);
			/*pot_en[j] = pot_en[j] - (G * masses[i] * masses[j])/(dpos[j][i][3]);*/
		}
	}
	return 0;
}


int Iterate(particle_i,timestep)										/*Particle iteration function*/
{
	data[timestep][particle_i][5] = pot_en[particle_i];					/*Add potetial energy to data array*/
	
	/*Calculate new velocity in each direction and add to new_vel array*/
	new_vel[particle_i][0] = vel[particle_i][0] + ((Fx[particle_i] * dt * 3600 * 24)/(masses[particle_i]));
	new_vel[particle_i][1] = vel[particle_i][1] + ((Fy[particle_i] * dt * 3600 * 24)/(masses[particle_i]));
	new_vel[particle_i][2] = vel[particle_i][2] + ((Fz[particle_i] * dt * 3600 * 24)/(masses[particle_i]));
	
	/*Calculate new x positon and add to data array*/
	new_pos[particle_i][0] = (1.0+a) * (pos[particle_i][0] + (new_vel[particle_i][0] * dt * 3600 * 24));
	/*printf ("Pos[i][0] - %.2e \t New_Pos[i][0] - %.2e \t New_Pos/(1.0+a) - %.2e \n",pos[particle_i][0],new_pos[particle_i][0],new_pos[particle_i][0]/(1.0+a));*/
	data[timestep][particle_i][0] = (new_pos[particle_i][0] / pow((1.0+a),(timestep+1)));
	
	/*Calculate new y positon and add to data array*/
	new_pos[particle_i][1] = (1.0+a) * (pos[particle_i][1] + (new_vel[particle_i][1] * dt * 3600 * 24));
	data[timestep][particle_i][1] = (new_pos[particle_i][1] / pow((1.0+a),(timestep+1)));
	
	/*Calculate new z positon and add to data array*/
	new_pos[particle_i][2] = (1.0+a) * (pos[particle_i][2] + (new_vel[particle_i][2] * dt * 3600 * 24));
	data[timestep][particle_i][2] = (new_pos[particle_i][2] / pow((1.0+a),(timestep+1)));
	
	/*Calculate velocities half a timestep for kinetic energy calcs*/
	double vel_x = vel[particle_i][0] + ((Fx[particle_i] * dt * 3600 * 24 * 0.5)/(masses[particle_i]));
	double vel_y = vel[particle_i][1] + ((Fy[particle_i] * dt * 3600 * 24 * 0.5)/(masses[particle_i]));
	double vel_z = vel[particle_i][2] + ((Fz[particle_i] * dt * 3600 * 24 * 0.5)/(masses[particle_i]));
	double vel_abs = sqrt(pow(vel_x,2) + pow(vel_y,2) + pow(vel_z,2));
	
	/*Calculate kinectic energy and add to kin_en array*/
	kin_en[particle_i] = 0.5*masses[particle_i]*pow(vel_abs,2);
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
		masses[particles] = M;
		pos[particles][0] = rand_range(min,max)*(size/2.0);
		pos[particles][1] = rand_range(min,max)*(size/2.0);
		pos[particles][2] = rand_range(min,max)*(size/2.0);
	}
	return 0;
}

int ImportParticles()
{
	FILE *particles_file;
	particles_file = fopen("particles.txt", "r");
	
	char line[100];
	int index = 0;
	
	fgets (line, sizeof line, particles_file); /*Skip first line*/
		
	while (fgets (line, sizeof line, particles_file))
	{
		//printf ("\n Line - %s \n", line);
		char *result = NULL;
		result = strtok(line, ",");
		masses[index] = atof(result);
		
		result = strtok( NULL, "," );
		pos[index][0] = atof(result);
		
		result = strtok( NULL, "," );
		pos[index][1] = atof(result);
		
		result = strtok( NULL, "," );
		pos[index][2] = atof(result);
		
		result = strtok( NULL, "," );
		vel[index][0] = atof(result);
		
		result = strtok( NULL, "," );
		vel[index][1] = atof(result);
		
		result = strtok( NULL, "," );
		vel[index][2] = atof(result);
		
		index++;
	}
	
	//printf ("\n x - %1.1e, y - %1.1e, z - %1.1e \t vx - %1.1e, vy - %1.1e, vz - %1.1e \n", pos[0][0],pos[0][1],pos[0][2],vel[0][0],vel[0][1],vel[0][2]);
	//printf ("\n x - %1.1e, y - %1.1e, z - %1.1e \t vx - %1.1e, vy - %1.1e, vz - %1.1e \n", pos[1][0],pos[1][1],pos[1][2],vel[1][0],vel[1][1],vel[1][2]);
	
	return 0;
}


/* #### Main Program #### */

int main()
{
	

	

	
	
	ImportParticles();
	//RandomPos(pos, masses, NumP, size, M); 														/*Asign random positions to particles*/
	
	/*
	pos[0][0] = 1.0E11;
	pos[1][0] = -1.0E11;
	pos[2][1] = 1.0E11;
	pos[3][1] = -1.0E11;
	
	vel[0][1] = 2.0E4;
	vel[1][1] = -2.0E4;
	vel[2][0] = -2.0E4;
	vel[3][0] = 2.0E4;

	masses[0] = 2E30;
	masses[1] = 2E30;
	masses[2] = 2E30;
	masses[3] = 2E30;
	*/
	
   	int t; 																/*Initiate timestep counter*/
   	for (t = 0; t < Ndt; t++) 											/*Loop through timesteps*/
   	{
		
		Calc();															/*Calculate dpos and force arrays for this timestep*/
		
    	int particle_i; 												/*Initiate particle counter*/
    	
    	for (particle_i = 0; particle_i < NumP; particle_i++) 			/*Loop through particles*/
    	{
			data[t][particle_i][3] = t * dt; 							/*Add time to data array*/
			Iterate(particle_i,t); 										/*Iterate particle_i*/
		}
		
		int i;
		
		for (i=0; i < NumP; i++) 										/*Iterate over particles*/
		{
			/*printf ("\n Posx %e \n", pos[i][0]);*/
			/*printf ("\n Posx/(1+a) %e \n", pos[i][0]/(1+a));*/
			int j;
			for (j=0; j < 3; j++)
			{
				pos[i][j] = new_pos[i][j]; 								/*Asign new positions to position array*/
				vel[i][j] = new_vel[i][j]; 								/*Asign new velocities to velocity array*/
			}
		}
   	}
   	
   	/*Data Output to stdio*/
   	
   	printf ("NumP,dt,Ndt,Precision, , ,\n");								
   	printf ("%i,%f,%i,%i,%e,%e,\n", NumP, dt, Ndt, prec, size, e);			/*System Parameters used in analysis*/
   	printf("x,y,z,t,Ek,Ep,\n");
   	int i;
   	
   	for (i=0; i < Ndt; i+=prec)											/*Loop through timesteps with specified precision*/
   	{
		
		int j;
		
		for (j=0; j < NumP; j++)										/*Loop over particles*/
		{
			/*Print particle data (x,y,z,t,Ek,Ep)*/
			printf("%e,%e,%e,%e,%e,%e,\n" , data[i][j][0], data[i][j][1], data[i][j][2], data[i][j][3], data[i][j][4], data[i][j][5]);
		}   	
   	}
   	   	   	
   	return 0;
}



