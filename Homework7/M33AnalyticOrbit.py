# Homework 7
# Analytic orbit of M33. Calculates the future orbit of M33 around M31 treating M33 as a point source, and using
# Hernquist and Miyamoto-Nagai acceleration profiles for M31.

# import modules
import numpy as np
import astropy.units as u
from astropy.constants import G
from ReadFile import Read
from CenterOfMass import CenterOfMass
import matplotlib
import matplotlib.pyplot as plt

class M33AnalyticOrbit:

    def __init__(self, filename):
        # initialize class, filename is the output file to contain the future locations and velocities of M33.
        # initial conditions are for M31 and M33 are loaded without input.

        #determine relative position and velocity of M33 com to M31 com, and save
        M31COM = CenterOfMass('/home/agibbs/M31_000.txt', 2)
        M33COM = CenterOfMass('/home/agibbs/M33_000.txt', 2)
        x31, y31, z31 = M31COM.COM_P(0.1)
        x33, y33, z33 = M33COM.COM_P(0.1)
        vx31, vy31, vz31 = M31COM.COM_V(x31, y31, z31)
        vx33, vy33, vz33 = M33COM.COM_V(x33, y33, z33)
        self.x, self.y, self.z = x33 - x31, y33 - y31, z33 - z31
        self.vx, self.vy, self.vz = vx33 - vx31, vy33 - vy31, vz33 - vz31

        # mass and characteristic radii for disk, bulge, and halo
        # to be used in Hernquist and Miyamoto-Nagai profiles
        self.rdisk = 5 * u.kpc
        self.Mdisk = 1.2e11 * u.Msun

        self.rbulge = 1 * u.kpc
        self.Mbulge = 1.9e10 * u.Msun

        self.rhalo = 60 * u.kpc
        self.Mhalo = 1.92e12 * u.Msun

        self.filename = filename

        return

    def HernquistAccel(self, m, ra, x, y, z, component):
        # returns one component of acceleration from a Hernquist profile given mass, characteristic radius,
        # relative positions, and acceleration component (x, y, or z)
        # valid for bulge and halo contributions
        if component == 'x':
            i = x
        elif component == 'y':
            i = y
        elif component == 'z':
            i=z
        else:
            print('Not a valid component')

        r = (x**2 + y**2 + z**2)**(0.5) # absolute distance of M33 to M31
        accel = -G * m / ( r * (ra + r)**2 ) * i
        return accel

    def MiyamotoNagaiAccel(self, m, rd, x, y, z, component):
        # returns one component of acceleration from a Miyamoto-Nagai profile given mass, characteristic radius,
        # relative positions, and acceleration component (x, y, or z)
        # valid for disk contributions
        zd = self.rdisk/5 #characteristic disk height
        B = rd + (z**2 + zd**2)**0.5
        R = (x**2 + y**2)**0.5
        accel = -G * m / (R**2 + B**2)**(1.5)
        if component == 'x':
            accel = accel * x
        elif component == 'y':
            accel = accel * y
        elif component == 'z': # z component is computed differently from x and y
            accel = accel * B * z / (z**2 + zd**2)**0.5
        else:
            print('Not a valid component')
        return accel

    def M31Accel(self, x, y, z, component):
        # returns one component of the summed acceleration from Hernquist (bulge, halo)
        # and Miyamoto-Nagai (disk) contributions given relative positions and component
        accel = self.MiyamotoNagaiAccel(self.Mdisk, self.rdisk, x, y, z, component) + \
                 self.HernquistAccel(self.Mbulge, self.rbulge, x, y, z, component) + \
                 self.HernquistAccel(self.Mhalo, self.rhalo, x, y, z, component)

        return accel

    def LeapFrog(self, dt, x, y, z, vx, vy, vz):
        # calculates future position and velocity of M33 a time dt from present, given current relative position and
        # velocity of M33

        # estimates position components of M33 at t + dt/2
        x1 = x + vx * dt / 2
        y1 = y + vy * dt / 2
        z1 = y + vz * dt / 2

        # estimates velocity components of M33 at t + dt, using positions at t + dt/2
        vx2 = vx + self.M31Accel(x1, y1, z1, 'x') * dt
        vy2 = vy + self.M31Accel(x1, y1, z1, 'y') * dt
        vz2 = vz + self.M31Accel(x1, y1, z1, 'z') * dt

        # estimates position components of M33 at t + dt by averaging velocity at t and t + dt
        x2 = x + 0.5 * (vx + vx2) * dt
        y2 = y + 0.5 * (vy + vy2) * dt
        z2 = z + 0.5 * (vz + vz2) * dt

        return x2, y2, z2, vx2, vy2, vz2

    def OrbitIntegrator(self, t0, dt, tmax):
        # calculates orbit of M33 from time t0 to tmax in steps of dt
        # by looping LeapFrog function
        # TIMES SHOULD BE INPUT IN YEARS

        # store initial conditions
        x = self.x
        y = self.y
        z = self.z
        vx = self.vx
        vy = self.vy
        vz = self.vz
        t = t0 * 31536000 * u.s #converts time in years to time in seconds
        dt = dt * 31536000  * u.s
        tmax = tmax * 31536000 * u.s

        # create initial arrays for each quantity, units are stripped for numpy
        xarr = np.array([x/u.kpc])
        yarr = np.array([y/u.kpc])
        zarr = np.array([z/u.kpc])
        vxarr = np.array([vx/u.km*u.s])
        vyarr = np.array([vy/u.km*u.s])
        vzarr = np.array([vz/u.km*u.s])
        tarr = np.array([t/u.s])

        while t < tmax:
            # iterate
            x, y, z, vx, vy, vz = self.LeapFrog(dt, x, y, z, vx, vy, vz)

            t = t + dt # update time

            # append arrays with new calculated quantities
            xarr = np.append(xarr, x/u.kpc)
            yarr = np.append(yarr, y/u.kpc)
            zarr = np.append(zarr, z/u.kpc)
            vxarr = np.append(vxarr, vx/u.km*u.s)
            vyarr = np.append(vyarr, vy/u.km*u.s)
            vzarr = np.append(vzarr, vz/u.km*u.s)
            tarr = np.append(tarr, t/u.s)

        # save orbit
        Orbit = np.array([tarr, xarr, yarr, zarr, vxarr, vyarr, vzarr])
        Orbit = Orbit.transpose()
        np.savetxt(self.filename, Orbit, header='t x y z vx vy vz', comments='#')

        return xarr, yarr, zarr, vxarr, vyarr, vzarr, tarr

# USING THE CLASS

#analytic orbit from this file
M33orbita = M33AnalyticOrbit('M33orbit.txt') #create class instance
x, y, z, vx, vy, vz, t = M33orbita.OrbitIntegrator(0.0, 0.1e9, 1e10)
r = (x**2 + y**2 + z**2)**0.5
v = (vx**2 + vy**2 + vz**2)**0.5
t = t / 3.1536e13 # seconds to Myr

#previous orbit from snapshots
dataM31 = np.genfromtxt('Orbit_M31.txt', dtype=None,names=True,skip_header=0)
dataM33= np.genfromtxt('Orbit_M33.txt', dtype=None,names=True,skip_header=0)
sep_M31_M33 = ( (dataM31['x'] - dataM33['x'])**2 + (dataM31['y'] - dataM33['y'])**2 + (dataM31['z'] - dataM33['z'])**2 )**0.5
vrel_M31_M33 = ( (dataM31['vx'] - dataM33['vx'])**2 + (dataM31['vy'] - dataM33['vy'])**2 + (dataM31['vz'] - dataM33['vz'])**2 )**0.5

time = dataM33['t']

#plotting
#plt.plot(t, r, 'r', time, sep_M31_M33, 'b')
plt.plot(t, v, 'r', time, vrel_M31_M33, 'b')
plt.xlabel('Time (Myr)')
plt.ylabel('Velocity (km/s)')
plt.title('Analytic Orbit of M33 around M31 - Velocity')
plt.show()
