# program adds MW and M31 files together
# used to create merger remnant galaxy files that can be used just like MW and M31 files

# import modules
import numpy as np
from ReadFile import Read

for i in range(0, 800, 10):
    # find MW and M31 files for snap
    ilbl = '000' + str(i) # add string of filenumber to 000
    ilbl = ilbl[-3:] # keep last 3 digits
    MWfile = '/home/agibbs/VLowRes/' + 'MW_' + ilbl +'.txt'
    M31file = '/home/agibbs/VLowRes/' + 'M31_' + ilbl +'.txt'
    # read file data
    MWt, MWtot, MWdata = Read(MWfile)
    M31t, M31tot, M31data = Read(M31file)
    # add file data together
    data = np.concatenate((MWdata, M31data), axis=0)
    t = 14.2857 * i
    fileout = '/home/agibbs/400B/ASTR400B_Gibbs/Project/Remnant/MW+M31_'+ilbl+'.txt'
    # save  remnant file with similar headers to original galaxy files
    np.savetxt(fileout, data, header='Time     %f \nTotal      135000 \n \
    mass in 1e10,  x, y, z, in kpc and vx, vy, vz in km/s \n type, m, x, y, z, vx, vy, vz'%t, comments='#')
