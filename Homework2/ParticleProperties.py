import numpy as np
import astropy.units as u
from ReadFile import Read


def ParticleInfo(particle_type, particle_num): #returns properties for any particle

    time, total, data = Read('/home/agibbs/MW_000.txt')

    #separate data types and particle types
    index = np.where(data['type'] == particle_type)

    x = np.around( data['x'][index] * u.kpc, 3 ) #only of particle type, keep 3 decimals, units kpc
    y = np.around( data['y'][index] * u.kpc, 3 )
    z = np.around( data['z'][index] * u.kpc, 3 )
    vx = np.around( data['vx'][index] * u.km / u.s, 3 )
    vy = np.around( data['vy'][index] * u.km / u.s, 3 )
    vz = np.around( data['vz'][index] * u.km / u.s, 3 )
    m = data['m'][index] * 1e10  * u.Msun

    n = particle_num #to make the next line shorter
    return x[n], y[n], z[n], vx[n], vy[n], vz[n], m[n]

x, y, z, vx, vy, vz, m = ParticleInfo(2, 99) #get 100th paticle of disk type
xn, yn, zn = x.to(u.lyr), y.to(u.lyr), z.to(u.lyr) #convert distance to light years

#print to terminal
print('x:', np.around(xn, 3), 'y:', np.around(yn, 3), 'z:', np.around(zn, 3)) #keep 3 decimals
print('vx:', vx, 'vy:', vy, 'vz:', vz)
print('Position magnitude:', np.around( (xn**2 + yn**2 + zn**2)**0.5, 3) )
print('Velocity magnitude:', np.around( (vx**2 + vy**2 + vz**2)**0.5, 3) )
print('mass:', m)
