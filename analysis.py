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
	Epsilon = int(data[1][5])
	
	print "Number of Particle -", NumP
	print "Timestep in Days -", dt
	print "Number of Timesteps -", Ndt
	print "Simulation Time -", "%0.1e" %(Ndt*dt/365.25), "years /", Ndt*dt, "days"
	print "The Universe size is %.1em" %Size
	print "The Epsilon Value is %.1em \n" %Epsilon
	
	return data,NumP,dt,Ndt,Prec,Size,Epsilon
	
def IterateData (Ndt,NumP,Prec,Size,Epsilon,plot):
	
	Ek_array = N.zeros(Ndt/Prec)
	Ep_array = N.zeros(Ndt/Prec)
	t_array = N.arange(0,Ndt/Prec,1)
	
	fig = plt.figure(num=None, figsize=(16, 9), dpi=80, facecolor='w', edgecolor='k')
	ax = p3.Axes3D(fig)
	
	t = 3
	t_actual = 0
	count = 0
	
	sys.stdout.write('Analysing Data...')
	
	while (t < (Ndt*NumP)/Prec):

		i = 0
		
		sys.stdout.write('.')
		sys.stdout.flush()
				
		area = float(N.pi * (Epsilon)**2)
				
		while (i < NumP):
			line = t+i
			
			if (plot == "2D"):
				plt.scatter(data[line][0],data[line][1],'bo',s=5)
			elif (plot == "3D"):
				ax.scatter(data[line][0],data[line][1],data[line][2],'bo',s=5)
									
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
	


def PlotEnergy(Ek,Ep,t,title,NumP,Ndt,dt):
	
	fig = plt.figure(num=None, figsize=(32, 18), dpi=60, facecolor='w', edgecolor='k')
	
	t = t * dt

	plt.plot(t,Ek,'r-', label='Kinetic')
	plt.plot(t,Ep,'g-', label='Potential')
	plt.plot(t,Ek+Ep,'b-', label='Total')
	plt.legend()
	
	E_init = Ek[0]+Ep[0]
	E_fin = Ek[-1]+Ep[-1]
	E_change = ((E_init - E_fin)/E_init)*100
	
	plttitle = 'No Particles: %1.1e' %(NumP) + ', Number of Timesteps: %1.1e' %(Ndt) + ', Timestep/days: %.3f' %(dt) + ', Initial Energy: %1.3e' %(E_init) + ' J, Final Energy: %1.3e' %(E_fin) + ' J, Percentage Change in Energy: %1.3e' %(E_change)
	name = 'Energy%s.png' %title
	plt.title(plttitle)
	plt.savefig(name)
	
def Video(fr, title):
	
	print "\n\nCreating Video ...\n"
	cmd = 'ffmpeg -r ' + str(fr) + ' -qscale 2 -i plots/Plot%05d.png ' + str(title) + '.mp4'
	subprocess.call(cmd, shell=True)
	

if __name__ == "__main__":

	parser = ap.ArgumentParser()
	parser.add_argument("title", help="Run ID", type=str)
	parser.add_argument("-d","--data", help="Location of csv data file. Default is 'data.txt.'", type=str)
	parser.add_argument("-p","--plots", help="Plot type - '2D', '3D'. If not called then no plots.", type=str)
	parser.add_argument("-e","--energy", help="Plot energy data. No argument required.", action="store_true")
	parser.add_argument("-v","--video", help="Generate Video Output. No argument required.", action="store_true")
	parser.add_argument("-f","--framerate", help="Video Framerate frames/s. Integer. Default 10fps", type=int)
	parser.add_argument("-n","--number", help="Number of Plots",type=int)
	args=parser.parse_args()
	
	if (args.title == None):
		print "No title/run ID defined. Use "
		sys.exit()
	else:
		print "\nRunning Analysis for Simulation - %s \n" %args.title
	
	if (args.data == None):
		data,NumP,dt,Ndt,Prec,Size,Epsilon = ImportData()
	else:
		data,NumP,dt,Ndt,Prec,Size,Epsilon = ImportData(args.data)

	Ek_array,Ep_array,t_array = IterateData(Ndt,NumP,Prec,Size,Epsilon,args.plots)
	
	PlotEnergy(Ek_array,Ep_array,t_array,args.title,NumP,Ndt,dt)
	
	if (args.plots != None and args.video != None):
		if (args.framerate == None):
			framerate = 10
		else:
			framerate = args.framerate
			
		Video(framerate,args.title)
	
	print "\nComplete\n"
	
	

