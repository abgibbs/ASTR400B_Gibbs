#Function to calculate the position of a galaxy from snap number start to end in increments of n
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from ReadFile import Read
from CenterOfMass import CenterOfMass


def OrbitCOM(galaxy, start, end, n):
    #setup paramters for CenterOfMass, array, and output file
    fileout = "Orbit_" + "%s"%(galaxy) +  ".txt"
    Orbit = np.zeros(( int(end/n) + 1, 7))
    delta = 0.5
    VolDec = 4

    #calculate center of mass for every snap number and save into Orbit array
    for i in np.arange(start, end+n, n):
        ilbl = '000' + str(i) # add string of filenumber to 000
        ilbl = ilbl[-3:] # keep last 3 digits
        filename = '/home/agibbs/VLowRes/' + "%s_"%(galaxy) + ilbl +'.txt'

        COM = CenterOfMass(filename, 2) #disk particles only
        xcom, ycom, zcom = COM.COM_P(delta, VolDec)
        vxcom, vycom, vzcom = COM.COM_V(xcom, ycom, zcom)

        Orbit[int(i/n), 0] = float(COM.time / u.Myr)
        Orbit[int(i/n), 1] = float(xcom / u.kpc)
        Orbit[int(i/n), 2] = float(ycom / u.kpc)
        Orbit[int(i/n), 3] = float(zcom / u.kpc)
        Orbit[int(i/n), 4] = float(vxcom / u.km * u.s)
        Orbit[int(i/n), 5] = float(vycom / u.km * u.s)
        Orbit[int(i/n), 6] = float(vzcom / u.km * u.s)

        print(i)

    #save text file with complete array
    np.savetxt(fileout, Orbit, header='t, x, y, z, vx, vy, vz', comments='#', fmt='%.2f')

    return

#OrbitCOM('M31', 0, 800, 5)

#PLOTTING
'''
#read data from all output files
dataMW = np.genfromtxt('Orbit_MW.txt', dtype=None,names=True,skip_header=0)
dataM31 = np.genfromtxt('Orbit_M31.txt', dtype=None,names=True,skip_header=0)
dataM33 = np.genfromtxt('Orbit_M33.txt', dtype=None,names=True,skip_header=0)

#calcualte relative positions and velocities
sep_MW_M31 = ( (dataMW['x'] - dataM31['x'])**2 + (dataMW['y'] - dataM31['y'])**2 + (dataMW['z'] - dataM31['z'])**2 )**0.5
sep_M31_M33 = ( (dataM31['x'] - dataM33['x'])**2 + (dataM31['y'] - dataM33['y'])**2 + (dataM31['z'] - dataM33['z'])**2 )**0.5
vrel_MW_M31 = ( (dataMW['vx'] - dataM31['vx'])**2 + (dataMW['vy'] - dataM31['vy'])**2 + (dataMW['vz'] - dataM31['vz'])**2 )**0.5
vrel_M31_M33 = ( (dataM31['vx'] - dataM33['vx'])**2 + (dataM31['vy'] - dataM33['vy'])**2 + (dataM31['vz'] - dataM33['vz'])**2 )**0.5

time = dataMW['t']

#plot
plt.plot(time, sep_M31_M33)
plt.title('Separation between M31 and M33')
plt.xlabel('Years from Present (Myr)')
plt.ylabel('Separation (kpc)')

plt.show()
'''
