
# Class to load halo and calculate spin and velocity anistropy
# Takes galaxy name as input and snaps to calculate kinematics for

# import modules
import numpy as np
import astropy.units as u
from astropy.constants import G
from ReadFile import Read
from CenterOfMass import CenterOfMass
#from r200 import r200
import matplotlib
import matplotlib.pyplot as plt

class HaloKinematics:

    def __init__(self, galaxy, start, end, n):
        # store snaps to calculate kinematics for
        self.galaxy = galaxy
        self.start = start
        self.n = n
        self.end = end

        # perform r200 calculations to match selected snaps, in case it has not been done
        #r200(galaxy, start, end, n)

        # load r200 values
        gal_r200 = np.genfromtxt('R200_'+galaxy+'.txt', dtype=None,names=True,skip_header=0)
        self.r200 = gal_r200['r200']
        self.qr200 = gal_r200['qr200']

        return

    def PosAndVels(self, filename, radius_cut):
        # calculates the relative positions and velocities of all halo particles, necessary for
        # vel anistropy and spin calculations
        # takes snap number and outer radius, returns position magnitude and components, and vel components

        # set COM params here
        delta = 5.0
        VolDec = 2

        COM = CenterOfMass(filename, 1)
        xcom, ycom, zcom = COM.COM_P(delta, VolDec)
        vxcom, vycom, vzcom = COM.COM_V(xcom, ycom, zcom)

        # relative positions
        x = COM.x - xcom
        y = COM.y - ycom
        z = COM.z - zcom

        R = ( x**2 + y**2 + z**2 )**0.5

        # relative velocities
        vx = COM.vx - vxcom
        vy = COM.vy - vycom
        vz = COM.vz - vzcom

        # only return particles within radius cut
        index = np.where(R < radius_cut * u.kpc)
        #print(radius_cut)

        return R[index], x[index], y[index], z[index], vx[index], vy[index], vz[index], COM.m[index] * u.Msun * 1e10


    def VelAnisotropy(self, filename, radius_cut):
        # calculates velocity anistropy given snap number and outer radius cut

        # get relative positions and velocities
        R, x, y, z, vx, vy, vz, m = self.PosAndVels(filename, radius_cut)

        # create radial unit vector
        xnorm = x / R
        ynorm = y / R
        znorm = z / R


        # dot product to find radial component squared
        vrad2 = (xnorm * vx + ynorm * vy + znorm * vz)**2

        # cross product to find tangential component squared
        vtan2 = (xnorm * vy - ynorm * vx)**2 + (-xnorm * vz + znorm * vx)**2 + (ynorm * vz - znorm *vy)**2

        vrsum = np.sum(vrad2)
        vtsum = np.sum(vtan2)
        vani = 1.0 - vtsum / vrsum / 2 # total velocity anistropy
        return vani

    def Spin(self, filename, radius_cut):
        # 'Practical' spin parameter from Bullocks 2001. spin = J/MVR(sqrt(2))
         R, x, y, z, vx, vy, vz, m = self.PosAndVels(filename, radius_cut)

         # calculate angular momentum vector
         Lz = np.sum((x * vy - y * vx) * m) * u.s / u.km / u.kpc / u.Msun
         Ly = np.sum((-x * vz + z * vx) * m) * u.s / u.km / u.kpc / u.Msun
         Lx = np.sum((y * vz - z *vy) * m) * u.s / u.km / u.kpc / u.Msun

         L = (Lz**2 + Ly**2 + Lx**2)**0.5

         Lxnorm = Lx / L
         Lynorm = Ly / L
         Lznorm = Lz / L

         # mass within radius cut
         M = np.sum(m)

         # circular velocity
         V = ( G * M / (radius_cut * u.kpc))**0.5 / u.km * u.s

         # spin
         spin = L / ((2**0.5) * M * V * radius_cut) * u.Msun
         #print(spin)
         return spin, Lxnorm, Lynorm, Lznorm

    def MergerProgression(self):
        # performs kinematic functions for each snap and saves values into text files

        # create arrays for data
        velani = np.zeros((int((self.end)/self.n) + 1, 3))
        spin200 = np.zeros((int((self.end)/self.n) + 1, 5))
        spinq200 = np.zeros((int((self.end)/self.n) + 1, 5))
        fileout1 = "VelANI_" + "%s"%(self.galaxy) +  ".txt" # velocity anisotropy file
        fileout2 = "Spin200_" + "%s"%(self.galaxy) +  ".txt" # spin file calculated out to r200
        fileout3 = "Spinq200_" + "%s"%(self.galaxy) +  ".txt" # spin file calculated out to 0.25 * r200

        counter = 0 # for r200 files

        for i in np.arange(self.start, self.end+self.n, self.n):
            # create filename for kinematic calculations
            ilbl = '000' + str(i) # add string of filenumber to 000
            ilbl = ilbl[-3:] # keep last 3 digits
            #filename = '/home/agibbs/VLowRes/' + "%s_"%(self.galaxy) + ilbl +'.txt' # change comments for MW+M31
            filename = '/home/agibbs/400B/ASTR400B_Gibbs/Project/Remnant/MW+M31_'+ilbl+'.txt'

            # perform velocity anisotropy calculations for this snap and add to array
            velani[int(i/self.n), 0] = i * 14.2857
            velani[int(i/self.n), 1] = self.VelAnisotropy(filename, self.r200[counter])
            velani[int(i/self.n), 2] = self.VelAnisotropy(filename, self.qr200[counter])

            # perform spin calculations and add to array
            L200, Lx200, Ly200, Lz200 = self.Spin(filename, self.r200[counter])
            Lq200, Lxq200, Lyq200, Lzq200 = self.Spin(filename, self.qr200[counter])
            print(L200)

            spin200[int(i/self.n), 0] = i * 14.2857
            spin200[int(i/self.n), 1] = L200
            spin200[int(i/self.n), 2] = Lx200
            spin200[int(i/self.n), 3] = Ly200
            spin200[int(i/self.n), 4] = Lz200

            spinq200[int(i/self.n), 0] = i * 14.2857
            spinq200[int(i/self.n), 1] = Lq200
            spinq200[int(i/self.n), 2] = Lxq200
            spinq200[int(i/self.n), 3] = Lyq200
            spinq200[int(i/self.n), 4] = Lzq200

            counter = counter + 1

        np.savetxt(fileout1, velani, header='t, r200_velANI, qr200_velANI', comments='#', fmt='%.4f')
        np.savetxt(fileout2, spin200, header='t, spin, Lx, Ly, Lz', comments='#', fmt='%.4f')
        np.savetxt(fileout3, spinq200, header='t, spin, Lx, Ly, Lz', comments='#', fmt='%.4f')

        return

# create class instance and calculate kinematic evolution
MWkin = HaloKinematics('MW+M31', 0, 800, 10)
MWkin.MergerProgression()
