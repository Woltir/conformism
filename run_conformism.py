import matplotlib.pyplot as plt
import numpy as np
from conformism import *
from time import time

seed = 1234
np.random.seed(seed)
T = np.linspace(0.1, 1., 64)
#T = [0.4]
V = 0.
#V = [0.1]
#A = 0.1/np.linspace(0.1, 15., 128)
#A = [.1, 0.5, 1.]
A = [1.]
N = 64
alpha= 0.1

tstep = 500

Energy = np.zeros(( len(T), len(A) ))
Magnetization = np.zeros(( len(T),len(A) ))
HeatCapacity = np.zeros(( len(T), len(A) ))
Susceptibility = np.zeros(( len(T), len(A) ))

tot = len(T)*len(A)
compt = 0

def total_time(t) :
        h, m, s = 0, 0, 0
        h = t/3600
        m = (t-3600*h)/60
        s = t - h*3600 - m*60
        return str(h) + ' h, ' + str(m) + ' min, ' + str(s) + ' s.'

print str(tot) + ' steps in total'

for i in range(len(T)) :
	for j in range(len(A)) :
		t1 = time()
		compt += 1
		c = conformism2D(N, T[i], A, V, alpha, tstep)
		c.simulation()
		t2 = time()
		if compt == 1:
			print str(t2-t1) + ' time for a step'
			print 'Expected simulation time : ' + total_time( int( (t2-t1)*tot))
		Energy[i,j] = sum(c.energy[-10:])/len(c.energy[-10:])
		Magnetization[i,j] = sum(c.magnetization[-10:])/len(c.magnetization[-10:])
		HeatCapacity[i,j] = sum(c.heatcapacity[-10:])/len(c.heatcapacity[-10:])
		Susceptibility[i,j] = sum(c.susceptibility[-10:])/len(c.susceptibility[-10:])
		print str(compt*100./tot) + ' % done'

f, axarr = plt.subplots(2, 2)
for k in range(len(A)) :
        axarr[0, 0].plot(T, Energy[:,k], label = 'A = ' + str(A[k]))
        axarr[0, 0].set_title('Energy vs Temperature')
        axarr[0, 1].plot(T, Magnetization[:,k], label = 'A = ' + str(A[k]))
        axarr[0, 1].set_title('Magnetization vs Temperature')
        axarr[1, 0].plot(T, HeatCapacity[:,k], label = 'A = ' + str(A[k]))
        axarr[1, 0].set_title('Heat Capacity vs Temperature')
        axarr[1, 1].plot(T, Susceptibility[:,k], label = 'A = ' + str(A[k]))
        axarr[1, 1].set_title('Susceptibility vs Temperature')
#plt.legend()
plt.show()


f2, ax = plt.subplots(1, 3)
ax[0].imshow(c.startconfig, interpolation = 'none')
ax[1].imshow(c.equilibriumconfig, interpolation = 'none')
ax[2].imshow(c.endconfig, interpolation = 'none')
plt.show()
