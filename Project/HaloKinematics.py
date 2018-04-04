# Homework 7
# Class to load halo and calculate spin and velocity anistropy
# Takes galaxy name as input and snaps to calculate kinematics for

# import modules
import numpy as np
import astropy.units as u
from astropy.constants import G
from ReadFile import Read
from CenterOfMass import CenterOfMass
from r200 import r200
import matplotlib
import matplotlib.pyplot as plt

class HaloKinematics:

    def __init__(self, galaxy, start, end, n):
        # store snaps to calculate kinematics for
        self.galaxy = galaxy
        self.start = start
        self.num = int(end/n)

        # perform r200 calculations to match selected snaps, in case it has not been done
        # r200(galaxy, start, end, n)

        # load r200 values
        gal_r200 = np.genfromtxt('R200_'+galaxy+'.txt', dtype=None,names=True,skip_header=0)
        self.r200 = gal_r200['r200']
        self.qr200 = gal_r200['qr200']

        return

    def PosAndVels(self, snap, radius_cut):
        # calculates the relative positions and velocities of all halo particles, necessary for
        # vel anistropy and spin calculations
        # takes snap number and outer radius, returns position magnitude and components, and vel components

        # set COM params here
        delta = 5.0
        VolDec = 2

        COM = CenterOfMass('/home/agibbs/VLowRes/' + "%s_"%(self.galaxy) + "%s"%(snap) +'.txt', 1)
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

        return R[index], x[index], y[index], z[index], vx[index], vy[index], vz[index]


    def VelAnistropy(self, snap, radius_cut):
        # calculates velocity anistropy given snap number and outer radius cut

        # get relative positions and velocities
        R, x, y, z, vx, vy, vz = self.PosAndVels(snap, radius_cut)

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

    def Spin(self):
        return

MWkin = HaloKinematics('MW', 0, 800, 20)
print(MWkin.VelAnistropy('000', MWkin.r200[0]))
