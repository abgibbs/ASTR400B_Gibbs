# Reads in angular momentum files for given galaxy and calculates angular change between vectors between each time step

# import modules
from ReadFile import Read
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

galaxy = 'MW'

galmom = np.genfromtxt('Spin200_'+galaxy+'.txt', dtype=None,names=True,skip_header=0)

t = galmom['t']
L = galmom['L']
Lx = galmom['Lx']
Ly = galmom['Ly']
Lz = galmom['Lz']


dotP = np.zeros(L.size-1)

for i in range(L.size-1):
    dotP[i] = Lx[i] * Lx[i+1] + Ly[i] * Ly[i+1] + Lz[i] * Lz[i+1]

dTheta = np.arccos(dotP) * 180 / np.pi
t = np.delete(t, 0, 0)
plt.plot(t, dTheta)
plt.title('Angular Change in Angular Momentum Vector')
plt.xlabel('Time (Myr)')
plt.ylabel('Change in Angle (Degrees)')
plt.show()
