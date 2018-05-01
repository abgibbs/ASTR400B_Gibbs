# plot velocity anisotropy evolution for MW, M31, and remnant from data files produced by HaloKinematics.py

import matplotlib.pyplot as plt
import numpy as np

#read data files
MW_velANI = np.genfromtxt('VelANI_MW.txt', dtype=None,names=True,skip_header=0)
M31_velANI = np.genfromtxt('VelANI_M31.txt', dtype=None,names=True,skip_header=0)
MWM31_velANI = np.genfromtxt('VelANI_MW+M31.txt', dtype=None,names=True,skip_header=0)
MW_pos = np.genfromtxt('Orbit_MW.txt', dtype=None,names=True,skip_header=0)
M31_pos = np.genfromtxt('Orbit_M31.txt', dtype=None,names=True,skip_header=0)

sep_MW_M31 = ( (MW_pos['x'] - M31_pos['x'])**2 + (MW_pos['y'] - M31_pos['y'])**2 + (MW_pos['z'] - M31_pos['z'])**2 )**0.5


# plotting 4 graphs, top is milky way va, M31 va, remnant va, bottom is galaxy sep
plt.figure(1)
# MW
plt.subplot(311)
plt.plot(MW_velANI['t'], MW_velANI['r200_velANI'], label='Within $r_{200}$')
plt.plot(MW_velANI['t'], MW_velANI['qr200_velANI'], 'r', label='Within $r_{200}/4$')
plt.ylabel('$\\beta$')
plt.axhline(y=0, color='black')
plt.axvline(x=4000, color='black', linestyle='--')
plt.axvline(x=5850, color='black', linestyle='--')
plt.axvline(x=6200, color='black', linestyle='-')
plt.legend(loc='upper right')
plt.title('Milky Way Velocity Anisotropy Evolution')
# M31
plt.subplot(312)
plt.plot(M31_velANI['t'], M31_velANI['r200_velANI'], label='Within $r_{200}$')
plt.plot(M31_velANI['t'], M31_velANI['qr200_velANI'], 'r', label='Within $r_{200}/4$')
plt.ylabel('$\\beta$')
plt.axhline(y=0, color='black')
plt.axvline(x=4000, color='black', linestyle='--')
plt.axvline(x=5850, color='black', linestyle='--')
plt.axvline(x=6200, color='black', linestyle='-')
plt.legend(loc='upper right')
plt.title('M31 Velocity Anisotropy Evolution')
# Remnant
plt.subplot(313)
plt.plot(MWM31_velANI['t'][40:], MWM31_velANI['r200_velANI'][40:], label='Within $r_{200}$') # only plot after merger
plt.plot(MWM31_velANI['t'][40:], MWM31_velANI['qr200_velANI'][40:], 'r', label='Within $r_{200}/4$')
plt.ylabel('$\\beta$')
plt.axhline(y=0, color='black')
plt.axvline(x=4000, color='black', linestyle='--')
plt.axvline(x=5850, color='black', linestyle='--')
plt.axvline(x=6200, color='black', linestyle='-')
plt.legend(loc='upper right')
plt.xlim(0, 12000)
plt.title('Remnant Velocity Anisotropy Evolution')
'''
# MW and M31 separation
plt.subplot(414)
plt.plot(MW_pos['t'], sep_MW_M31)
plt.xlabel('Time (Myr)')
plt.ylabel('MW-M31 sep (kpc)')
plt.title('MW-M31 Distance')
'''
plt.show()
