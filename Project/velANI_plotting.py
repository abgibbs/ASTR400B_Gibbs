import matplotlib.pyplot as plt
import numpy as np

#read data files
MW_velANI = np.genfromtxt('VelANI_MW.txt', dtype=None,names=True,skip_header=0)
M31_velANI = np.genfromtxt('VelANI_M31.txt', dtype=None,names=True,skip_header=0)
MW_pos = np.genfromtxt('Orbit_MW.txt', dtype=None,names=True,skip_header=0)
M31_pos = np.genfromtxt('Orbit_M31.txt', dtype=None,names=True,skip_header=0)

sep_MW_M31 = ( (MW_pos['x'] - M31_pos['x'])**2 + (MW_pos['y'] - M31_pos['y'])**2 + (MW_pos['z'] - M31_pos['z'])**2 )**0.5


# plotting 3 graphs, top is milky way va, mid is M31 va, bottom is galaxy sep
plt.figure(1)
plt.subplot(311)
plt.plot(MW_velANI['t'], MW_velANI['r200_velANI'], label='Within r200')
plt.plot(MW_velANI['t'], MW_velANI['qr200_velANI'], 'r', label='Within 0.25xr200')
plt.ylabel('Vel Anisotropy')
plt.axhline(y=0, color='black')
plt.legend(loc='upper right')
plt.title('Milky Way Velocity Anistropy Evolution')

plt.subplot(312)
plt.plot(M31_velANI['t'], M31_velANI['r200_velANI'], label='Within r200')
plt.plot(M31_velANI['t'], M31_velANI['qr200_velANI'], 'r', label='Within 0.25xr200')
plt.ylabel('Vel Anisotropy')
plt.axhline(y=0, color='black')
plt.legend(loc='upper right')
plt.title('M31 Velocity Anistropy Evolution')

plt.subplot(313)
plt.plot(MW_pos['t'], sep_MW_M31)
plt.xlabel('Time (Myr)')
plt.ylabel('MW-M31 sep (kpc)')
plt.title('MW-M31 Distance')

plt.show()
