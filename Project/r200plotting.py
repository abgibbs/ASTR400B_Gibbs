# plot the r200 evolution of the MW, M31, and remnant as calculated and saved by r200.py

import matplotlib.pyplot as plt
import numpy as np

#read data files
MW_r200 = np.genfromtxt('R200_MW.txt', dtype=None,names=True,skip_header=0)
M31_r200 = np.genfromtxt('R200_M31.txt', dtype=None,names=True,skip_header=0)
MWM31_r200 = np.genfromtxt('R200_MW+M31.txt', dtype=None,names=True,skip_header=0)
MW_pos = np.genfromtxt('Orbit_MW.txt', dtype=None,names=True,skip_header=0)
M31_pos = np.genfromtxt('Orbit_M31.txt', dtype=None,names=True,skip_header=0)

sep_MW_M31 = ( (MW_pos['x'] - M31_pos['x'])**2 + (MW_pos['y'] - M31_pos['y'])**2 + (MW_pos['z'] - M31_pos['z'])**2 )**0.5

#plotting 2 graphs, top is r200 halo radius, bottom is galaxy sep

plt.figure(1)
# r200
plt.subplot(211)
plt.plot(MW_r200['t'], MW_r200['r200'], label='MW')
plt.ylabel(' $r_{200}$ (kpc)')
plt.title('$r_{200}$ Evolution of MW, M31, and Remnant')

#plt.subplot(412)
plt.plot(M31_r200['t'], M31_r200['r200'], 'r', label='M31')
#plt.ylabel('M31 r200 (kpc)')

#plt.subplot(413)
plt.plot(MWM31_r200['t'][40:], MWM31_r200['r200'][40:], 'g', label='Remnant')
#plt.ylabel('MW+M31 Remnant r200 (kpc)')
plt.legend(loc='upper left')
# MW and M31 separation
plt.subplot(212)
plt.plot(MW_pos['t'], sep_MW_M31)
plt.xlabel('Time (Myr)')
plt.ylabel('MW-M31 sep (kpc)')

plt.show()
