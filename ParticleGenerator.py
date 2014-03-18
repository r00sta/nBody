### Particle Generator Program foo nBody ###

import numpy as N	#Maths modules
import sys	#System controls
import random

def GlobCluster (numP, cenx, ceny, cenz, posrange, velrange, massmin, massmax, particle_file) :
	print "Making a Globular Cluster"
	
	i = 0
	
	while (i < numP):
		
		mass = 1.0E0 + random.randrange(massmin,massmax,massmin/1000000)
		posx = cenx + random.randrange(-posrange,posrange,posrange/1000000)
		posy = ceny + random.randrange(-posrange,posrange,posrange/1000000)
		posz = cenz + random.randrange(-posrange,posrange,posrange/1000000)
		velx = 0.0#random.randrange(-velrange,velrange,velrange/10)
		vely = 1.0e+7#random.randrange(-velrange,velrange,velrange/10)
		velz = 0.0#random.randrange(-velrange,velrange,velrange/10)
		
		line = str(mass) + "," + str(posx) + "," + str(posy) + "," + str(posz) + "," + str(velx) + "," + str(vely) + "," + str(velz) + "\n"
		
		print line
		
		particle_file.write(line)
		
		i+=1
	
	



if __name__ == "__main__":
	
	print "\nParticle Generator\n"

	particles_file = open("particles.txt", "w")
	
	particles_file.write("Mass,Posx,Posy,Posz,Velx,Vely,Velz\n")

	particles_file.write("1.0E40,0.0,0.0,0.0,0.0,0.0,0.0\n")
	
	GlobCluster(499,5.0E15,0.0,0.0,1.0E15,1.0E2,1.0E27,1.0E30,particles_file)
	
	#GlobCluster(100,5.0E17,0.0,0.0,1.0E17,1.0E3,1.0E27,1.0E34,particles_file)

	particles_file.close()
