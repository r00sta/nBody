### Data Analysis Program for nBody simulator. ###

import numpy as N	#Maths modules
import sys	#System controls
import subprocess
import matplotlib.pyplot as plt #Matplotlib module
import matplotlib.animation as animation
import mpl_toolkits.mplot3d.axes3d as p3	#3D axis
import argparse as ap	#Argument parser

def ImportData (datafile="data.txt"):
	
	print "Reading in data file %s ...\n" %datafile
	
	data = N.genfromtxt(datafile, delimiter=',')
	
	NumP = int(data[1][0])
	dt = float(data[1][1])
	Ndt = int(data[1][2])
	Prec = int(data[1][3])
	Size = int(data[1][4])
	
	print "Number of Particle -", NumP
	print "Timestep in Days -", dt
	print "Number of Timesteps -", Ndt
	print "Simulation Time -", "%0.1e" %(Ndt*dt/365.25), "years /", Ndt*dt, "days"
	print "The Universe size is %.1em \n" %Size
	
	return data,NumP,dt,Ndt,Prec,Size
	
def Iterate(Ndt,NumP,Prec,Size,plot):
	
	Ek_array = N.zeros(Ndt/Prec)
	Ep_array = N.zeros(Ndt/Prec)
	t_array = N.arange(0,Ndt/Prec,1)
	
	t = 3
	t_actual = 0
	count = 0
	
	sys.stdout.write('Analysing Data...')
	
	while (t < (Ndt*NumP)/Prec):

		i = 0
		
		sys.stdout.write('.')
		sys.stdout.flush()
				
		while (i < NumP):
			line = t+i
			
			if (plot == "2D"):
				plt.scatter(data[line][0],data[line][1],'bo')
			elif (plot == "3D"):
				ax.scatter(data[line][0],data[line][1],data[line][2],'bo')
						
			Ek_array[t_actual] += data[line][4]
			Ep_array[t_actual] += data[line][5]
		
			i+=1
		
		if (plot == "2D"):
			plt.xlim(-Size,Size)
			plt.ylim(-Size,Size)
			
		elif (plot == "3D"):
			ax.set_xlim3d([-Size, Size])
			ax.set_xlabel('X')
			ax.set_ylim3d([-Size, Size])
			ax.set_xlabel('Y')
			ax.set_zlim3d([-Size, Size])
			ax.set_xlabel('Z')
			
		#ax.set_frame_on(False)
		#ax.grid(False)
		#ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
		#ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
		#ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
		#ax.grid(True)
		#ax.set_axis_off()
		
		name = 'plots/Plot' + "%05d" %(count) + '.png'
		count+=1
		time = 'Time - ' + str(data[line][3]) + ' days '
		plt.title(time)
		plt.savefig(name)
		plt.cla()
		t+=(NumP)
		t_actual += 1
		
	return Ek_array,Ep_array,t_array
	
if __name__ == "__main__":

	data,NumP,dt,Ndt,Prec,Size = ImportData()

	Ek_array,Ep_array,t_array = IterateData(Ndt,NumP,Prec,Size,args.plots)

	fig = plt.figure()
	
	ani = animation.FuncAnimation(fig, Iterate, interval=50, blit=True)
	
	plt.show()
