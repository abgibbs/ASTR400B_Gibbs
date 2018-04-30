import numpy as np
import astropy.units as u

def Read(filename): #read header and save data from file
    file = open(filename, 'r')

    #read the time in Myr from the first line of the file
    line1 = file.readline()
    label1, value1 = line1.split()
    time = float(value1) * u.Myr

    #read the total particles from the second line of file
    line2 = file.readline()
    label2, value2 = line2.split()
    total = float(value2)

    file.close()

    #save the data from file for future access
    data = np.genfromtxt(filename, dtype=None,names=True,skip_header=3)

    return time, total, data

