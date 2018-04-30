# plot spin evolution of MW, M31, and remnant using data files produced by HaloKinematics.py

import matplotlib.pyplot as plt
import numpy as np

#read data files
MW_s200 = np.genfromtxt('Spin200_MW.txt', dtype=None,names=True,skip_header=0)
MW_qs200 = np.genfromtxt('Spinq200_MW.txt', dtype=None,names=True,skip_header=0)
M31_s200 = np.genfromtxt('Spin200_M31.txt', dtype=None,names=True,skip_header=0)
M31_qs200 = np.genfromtxt('Spinq200_M31.txt', dtype=None,names=True,skip_header=0)
MWM31_s200 = np.genfromtxt('Spin200_MW+M31.txt', dtype=None,names=True,skip_header=0)
MWM31_qs200 = np.genfromtxt('Spinq200_MW+M31.txt', dtype=None,names=True,skip_header=0)
MW_pos = np.genfromtxt('Orbit_MW.txt', dtype=None,names=True,skip_header=0)
M31_pos = np.genfromtxt('Orbit_M31.txt', dtype=None,names=True,skip_header=0)

sep_MW_M31 = ( (MW_pos['x'] - M31_pos['x'])**2 + (MW_pos['y'] - M31_pos['y'])**2 + (MW_pos['z'] - M31_pos['z'])**2 )**0.5

# plot 4 graphs, top MW, M31, remnant, bottom is MW and M31 separation

plt.figure(1)
# MW
plt.subplot(411)
plt.plot(MW_s200['t'], MW_s200['spin'], label='Within $r_{200}$')
plt.plot(MW_qs200['t'], MW_qs200['spin'], 'r', label='Within $r_{200}/4$')
plt.ylabel("Spin $\lambda'$")
plt.legend(loc='upper left')
plt.title('MW Spin Evolution - Within $r_{200}$')
# M31
plt.subplot(412)
plt.plot(M31_s200['t'], M31_s200['spin'], label='Within $r_{200}$')
plt.plot(M31_qs200['t'], M31_qs200['spin'], 'r', label='Within $r_{200}/4$')
plt.ylabel("Spin $\lambda'$")
plt.legend(loc='upper left')
plt.title('M31 Spin Evolution - Within $r_{200}$')
# Remnant
plt.subplot(413)
plt.plot(MWM31_s200['t'][40:], MWM31_s200['spin'][40:], label='Within $r_{200}$') #only plot after galaxies merge
plt.plot(MWM31_qs200['t'][40:], MWM31_qs200['spin'][40:], 'r', label='Within $r_{200}/4$')
plt.xlim(0, 12000)
plt.ylabel("Spin $\lambda'$")
plt.legend(loc='upper left')
plt.title('Remnant Spin Evolution - Within $r_{200}$')
# MW and M31 separation
plt.subplot(414)
plt.plot(MW_pos['t'], sep_MW_M31)
plt.xlabel('Time (Myr)')
plt.ylabel('MW-M31 sep (kpc)')
plt.title('MW-M31 Distance')

plt.show()
