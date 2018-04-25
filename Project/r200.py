# Function to calculate the r200, and 0.25*r200 radius given galaxy, and snap numbers from start to end in n spaces
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from ReadFile import Read
from CenterOfMass import CenterOfMass


def r200(galaxy, start, end, n):
    #setup paramters for CenterOfMass, array, and output file
    fileout = "R200_" + "%s"%(galaxy) +  ".txt"
    R200 = np.zeros(( int(end/n) + 1, 3))
    delta = 5.0
    VolDec = 2
    crit_density = 147.712 * u.Msun / u.kpc**3

    #calculate r200 for every snap number and save into R200 array
    for i in np.arange(start, end+n, n):
        ilbl = '000' + str(i) # add string of filenumber to 000
        ilbl = ilbl[-3:] # keep last 3 digits
        filename = '/home/agibbs/VLowRes/' + "%s_"%(galaxy) + ilbl +'.txt'

        #COM position
        COM = CenterOfMass(filename, 1) #halo particles only
        xcom, ycom, zcom = COM.COM_P(delta, VolDec)

        #Find R200 radius
        xref = COM.x - xcom
        yref = COM.y - ycom
        zref = COM.z - zcom
        R = (xref**2 + yref**2 + zref**2)**0.5 
        density = 1e9 * u.Msun / u.kpc**3 # initial density to start loop
        R_step = 0.1 * u.kpc # radius steps
        R_iter = 50.0 * u.kpc # first radius
        while density > crit_density * 200:
            R_iter = R_iter + R_step
            rindex = np.where(R < R_iter)
            mass = np.sum(COM.m[rindex]) * 1e10 * u.Msun
            volume = 4*np.pi/3 * R_iter**3
            density = mass / volume



        #Save
        R200[int(i/n), 0] = float(COM.time / u.Myr)
        R200[int(i/n), 1] = float(R_iter / u.kpc)
        R200[int(i/n), 2] = float(R_iter * 0.25 / u.kpc)

        print(i)

    #save text file with complete array
    np.savetxt(fileout, R200, header='t, r200, qr200', comments='#', fmt='%.2f')

    return

r200('MW', 0, 800, 10)
