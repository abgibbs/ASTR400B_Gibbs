from ReadFile import Read
import astropy.units as u
import numpy as np

#computes total mass of particle type in galaxy given by data in filename

def ComponentMass(filename, particle_type):

    time, total, data =  Read(filename)
    index = np.where(data['type'] == particle_type) #select only particle type
    m_tot = np.sum(data['m'][index]) * 1e10 * u.Msun #sum all masses of particle type, put into units of solar mass

    return np.around(m_tot, 3) #round to 3 decimal places

print(ComponentMass('/home/agibbs/M33_000.txt', 3))
