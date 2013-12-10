### Data Analysis Program for n_body2 simulator. ###
import numpy as N
import scipy as S
import pylab as py
import random as R
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3

data = N.genfromtxt('data.txt', delimiter=',')

NumP = int(data[1][0])
dt = float(data[1][1])
Ndt = int(data[1][2])
Prec = int(data[1][3])
Size = int(data[1][4])

print "\nNumber of Particle -", NumP
print "Timestep in Days -", dt
print "Number of Timesteps -", Ndt
print "Simulation Time -", "%0.1e" %(Ndt*dt/365.25), "years /", Ndt*dt, "days\n"

t = 3
count = 0
t_actual = 0
scale = 1


Ek_array = N.zeros(Ndt/Prec)
Ep_array = N.zeros(Ndt/Prec)
t_array = N.arange(0,Ndt/Prec,1)


#py.figure(num=None, figsize=(16, 9), dpi=80, facecolor='w', edgecolor='k')

fig = plt.figure(num=None, figsize=(16, 9), dpi=80, facecolor='w', edgecolor='k')
ax = p3.Axes3D(fig)



while (t < (Ndt*NumP)/Prec):
	#print t_actual
	i = 0
	print t
	#t_array[t_actual] = t/NumP
	
	while (i < NumP):
		line = t+i
		#print data[line]
		### Position
		#print line
		
		#py.plot(data[line][0],data[line][1],'bo')
		ax.scatter(data[line][0],data[line][1],data[line][2],'bo')
		
		
		Ek_array[t_actual] += data[line][4]
		Ep_array[t_actual] += data[line][5]
		
		i+=1
	
	#print t_array[t_actual],Ek_array[t_actual],t_actual
	
	#py.xlim(-2E11,2E11)
	#py.ylim(-2E11,2E11)
	
	ax.set_xlim3d([-Size, Size])
	ax.set_xlabel('X')
	ax.set_ylim3d([-Size, Size])
	ax.set_xlabel('Y')
	ax.set_zlim3d([-Size, Size])
	ax.set_xlabel('Z')
	
	name = 'plots/Plot' + "%05d" %(count) + '.png'
	count+=1
	time = 'Time - ' + str(data[line][3]) + ' days '
	plt.title(time)
	plt.savefig(name)
	plt.cla()
	
	#py.show()
	t+=(NumP)
	t_actual += 1
	
### Energy Plot

#print t_array

#print t_array[-2],t_array[-1]

py.figure(num=None, figsize=(32, 18), dpi=60, facecolor='w', edgecolor='k')

t_array = t_array * dt / scale

py.plot(t_array,Ek_array,'r-', label='Kinetic')
py.plot(t_array,Ep_array,'g-', label='Potential')
py.plot(t_array,Ek_array+Ep_array,'b-', label='Total')
py.legend()
#title = 'Energy'
E_init = Ek_array[0]+Ep_array[0]
E_fin = Ek_array[-1]+Ep_array[-1]
E_change = ((E_init - E_fin)/E_init)*100

title = 'No Particles: %1.1e' %(NumP) + ', Number of Timesteps: %1.1e' %(Ndt) + ', Timestep/days: %.3f' %(dt) + ', Initial Energy: %1.3e' %(E_init) + ' J, Final Energy: %1.3e' %(E_fin) + ' J, Percentage Change in Energy: %1.3e' %(E_change)
name = 'Energy.png'
py.title(title)
py.savefig(name)
#py.xlim(0,Ndt*dt)
#py.ylim(-4E39,4E39)
#py.show()
	
