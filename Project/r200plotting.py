import matplotlib.pyplot as plt
import numpy as np

#read data files
MW_r200 = np.genfromtxt('R200_MW.txt', dtype=None,names=True,skip_header=0)
M31_r200 = np.genfromtxt('R200_M31.txt', dtype=None,names=True,skip_header=0)
MW_pos = np.genfromtxt('Orbit_MW.txt', dtype=None,names=True,skip_header=0)
M31_pos = np.genfromtxt('Orbit_M31.txt', dtype=None,names=True,skip_header=0)

sep_MW_M31 = ( (MW_pos['x'] - M31_pos['x'])**2 + (MW_pos['y'] - M31_pos['y'])**2 + (MW_pos['z'] - M31_pos['z'])**2 )**0.5

#plotting 3 graphs, top is MW halo radius, mid is M31 halo radius, bottom is galaxy sep

plt.figure(1)
plt.subplot(311)
plt.plot(MW_r200['t'], MW_r200['r200'])
plt.ylabel('MW r200 (kpc)')

plt.subplot(312)
plt.plot(M31_r200['t'], M31_r200['r200'])
plt.ylabel('M31 r200 (kpc)')

plt.subplot(313)
plt.plot(MW_pos['t'], sep_MW_M31)
plt.xlabel('Time (Myr)')
plt.ylabel('MW-M31 sep (kpc)')

plt.show()
