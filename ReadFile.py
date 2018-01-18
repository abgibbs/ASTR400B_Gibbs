import numpy as np
import astropy.units as u

def Read(filename):
    file = open(filename, 'r')
    line1 = file.readline()
    label1, value1 = line1.split()
    time = float(value1) * 10.0 * u.Myr

    line2 = file.readline()
    label2, value2 = line2.split()
    total = float(value2)

    file.close()

    data = np.genfromtxt(filename, dtype=None,names=True,skip_header=3)
    
    return time, total, data



