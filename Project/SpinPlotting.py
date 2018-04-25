import matplotlib.pyplot as plt
import numpy as np

#read data files
gal_s200 = np.genfromtxt('Spin200_M31.txt', dtype=None,names=True,skip_header=0)
gal_qs200 = np.genfromtxt('Spinq200_M31.txt', dtype=None,names=True,skip_header=0)
MW_pos = np.genfromtxt('Orbit_MW.txt', dtype=None,names=True,skip_header=0)
M31_pos = np.genfromtxt('Orbit_M31.txt', dtype=None,names=True,skip_header=0)

sep_MW_M31 = ( (MW_pos['x'] - M31_pos['x'])**2 + (MW_pos['y'] - M31_pos['y'])**2 + (MW_pos['z'] - M31_pos['z'])**2 )**0.5


plt.figure(1)
plt.subplot(311)
plt.plot(gal_s200['t'], gal_s200['L'], label='Within r200')
plt.ylabel('Ang. Mom. (Msun kpc km/s)')
plt.title('M31 Angular Momentum Evolution - Within r200')

plt.subplot(312)
plt.plot(gal_qs200['t'], gal_qs200['L'], 'r', label='Within 0.25xr200')
plt.ylabel('Ang. Mom. (Msun kpc km/s)')
plt.title('Within 0.25 x r200')

plt.subplot(313)
plt.plot(MW_pos['t'], sep_MW_M31)
plt.xlabel('Time (Myr)')
plt.ylabel('MW-M31 sep (kpc)')
plt.title('MW-M31 Distance')

plt.show()
