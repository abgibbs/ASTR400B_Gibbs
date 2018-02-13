# Homework 5
# Class to determine mass and velocity for different particle types for any given galaxy and time.

# import modules
import numpy as np
import astropy.units as u
from astropy.constants import G
from ReadFile import Read
from CenterOfMass import CenterOfMass
import matplotlib
import matplotlib.pyplot as plt


class MassProfile:

    def __init__(self, galaxy, snap):
        # construct filename
        ilbl = '000' + str(snap) # add string of filenumber to 000
        ilbl = ilbl[-3:] # keep last 3 digits
        self.filename = '/home/agibbs/' + "%s_"%(galaxy) + ilbl +'.txt'

        # set delta for all com calculations
        self.comdel = 0.3 # too small a value may fail to converge and can result in an ERROR!!!

        # set G to correct units
        self.G = G.to(u.kpc*u.km**2/u.s**2/u.Msun)

        # store galaxy name (used to separate out M33 later)
        self.gname = galaxy

        # read in the file and particle type
        self.time, self.total, self.data = Read(self.filename)

        # store the mass, positions of all particle types
        self.m = self.data['m']

        self.x = self.data['x'] * u.kpc
        self.y = self.data['y'] * u.kpc
        self.z = self.data['z'] * u.kpc

    def MassEnclosed(self, ptype, radii):
        # Return array of masses within specific radii of COM. Takes particle type and an array of radii to use.
        p_index = np.where(self.data['type'] == ptype)
        GALCOM = CenterOfMass(self.filename, ptype)
        GAL_xcom, GAL_ycom, GAL_zcom = GALCOM.COM_P(self.comdel)
        radius = ( (self.x[p_index] - GAL_xcom)**2 + (self.y[p_index] - GAL_ycom)**2 + \
                   (self.z[p_index] - GAL_zcom)**2 ) ** 0.5
        mass = np.zeros(radii.size)
        p_mass = self.m[p_index]
        #print(np.argwhere(np.isnan(radius))) #debugging nan values
        for index, radii_val in enumerate(radii):
            r_index = np.where(radius < radii_val)
            mass[index] = np.sum(p_mass[r_index])
        return mass * u.Msun * 1e10

    def TotalMassEnclosed(self, radii):
        # Performs MassEnclosed for all particles added together. Takes array of radii to use.
        if self.gname != 'M33': #M33 does not have a bulge
            mass3 = self.MassEnclosed(3, radii)
        mass2 = self.MassEnclosed(2, radii)
        mass1 = self.MassEnclosed(1, radii)
        try:
            mass3
        except NameError:
            m_tot = mass2 + mass1
        else:
            m_tot = mass1 + mass2 + mass3
        return m_tot

    def HernquistMass(self, radii, a, Mhalo):
        # Calculates an array of mass enclosed using a theoretical model given an array of to use,
        # scale factor a, and halo mass.
        hmass = Mhalo * radii**2 / (a + radii)**2
        return hmass

    def CircularVelocity(self, ptype, radii):
        # Calculates the circular velocity for particle by balancing gravity from enclosed mass and centripetal force.
        # Takes array of radii to use and particle type.
        masses = self.MassEnclosed(ptype, radii)
        velocity = ( self.G * masses / radii)**0.5
        return np.around(velocity, 2)

    def TotalCircularVelocity(self, radii):
        # Performs circular velocity for all particles, using total mass enclosed.
        # Takes array of radii to use.
        masses = self.TotalMassEnclosed(radii)
        velocity = ( self.G * masses / radii )**0.5
        return np.around(velocity, 2)

    def HernquistVCirc(self, radii, a, Mhalo):
        # Calculates theoretical circular velocity using herniquist mass enclosed.
        # Takes hernquist parameters and array of radii to use.
        masses = self.HernquistMass(radii, a, Mhalo)
        velocity = ( self.G * masses / radii )**0.5
        return np.around(velocity, 2)

# CALCULATIONS
GALPROF = MassProfile('M33', 000)
radii = np.linspace(0.1, 100, 400) * u.kpc
a = 15.0 * u.kpc

# Mass Profile
masses1 = GALPROF.MassEnclosed(1, radii)
masses2 = GALPROF.MassEnclosed(2, radii)
#masses3 = GALPROF.MassEnclosed(3, radii)
m_tot = GALPROF.TotalMassEnclosed(radii)
hmass = GALPROF.HernquistMass(radii, a, masses1[-1])

# Rotation Curve
vel1 = GALPROF.CircularVelocity(1, radii)
vel2 = GALPROF.CircularVelocity(2, radii)
#vel3 = GALPROF.CircularVelocity(3, radii)
v_tot = GALPROF.TotalCircularVelocity(radii)
hvel = GALPROF.HernquistVCirc(radii, a, masses1[-1])

# PLOTTING

fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)

# Mass profile
'''
plt.semilogy(radii, masses1, color='blue', label='Dark Matter')
plt.semilogy(radii, masses2, color='red', label='Disk')
#plt.semilogy(radii, masses3, color='green', label='Bulge')
plt.semilogy(radii, m_tot, color='black', label='Total')
plt.semilogy(radii, hmass, color='cyan', label='Hernquist (a=15.0 kpc)', linestyle='--')

plt.xlabel('Radius (kpc)')
plt.ylabel('Solar Masses')

plt.xlim(0, 30)

plt.title('M33 Mass Profile')
'''
# Rotation Curve

plt.plot(radii, vel1, color='blue', label='Dark Matter')
plt.plot(radii, vel2, color='red', label='Disk')
#plt.plot(radii, vel3, color='green', label='Bulge')
plt.plot(radii, v_tot, color='black', label='Total')
plt.plot(radii, hvel, color='cyan', label='Hernquist (a=15.0 kpc)', linestyle='--')

plt.xlabel('Radius (kpc)')
plt.ylabel('Velocity (km/s)')

plt.xlim(0, 30)

plt.title('M33 Rotation Curve')

legend = ax.legend(loc='lower right')

plt.show()
