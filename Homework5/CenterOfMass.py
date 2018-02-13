# Homework 4
# Center of Mass Position and Velocity

# import modules
import numpy as np
import astropy.units as u
from ReadFile import Read


class CenterOfMass:

    def __init__(self, filename, ptype):
        # read in the file and particle type
        self.time, self.total, self.data = Read(filename)

        #create an array to store indexes of particles of desired Ptype
        self.index = np.where(self.data['type'] == ptype)

        # store the mass, positions, velocities of only the particles of the given type
        self.m = self.data['m'][self.index]

        ##### PLACE other particle properties here: x,y,z,vx,vy,vz #####
        self.x = self.data['x'][self.index] * u.kpc
        self.y = self.data['y'][self.index] * u.kpc
        self.z = self.data['z'][self.index] * u.kpc

        self.vx = self.data['vx'][self.index] * u.km / u.s
        self.vy = self.data['vy'][self.index] * u.km / u.s
        self.vz = self.data['vz'][self.index] * u.km / u.s

    def total_mass(self):
        #Note: you can add other keyword arguments into the function, but 'self' must be first
        return np.sum(self.m)*u.Msun*1e10

        ##### PLACE OTHER FUNCTIONS BELOW #####

    def COMdefine(self, x, y, z, index):
        # x, y, z refers to coordinate but could be position or velocity, index chooses which particles to use
        x_com = np.sum(self.m[index] * x[index]) / np.sum(self.m[index]) # definition of COM
        y_com = np.sum(self.m[index] * y[index]) / np.sum(self.m[index])
        z_com = np.sum(self.m[index] * z[index]) / np.sum(self.m[index])
        return x_com, y_com, z_com

    def COM_P(self, delta):
        # Compute center of mass position first with all particles, then with smaller volumes until convergence
        # Delta sets maximum COM change required for convergence
        XCOM, YCOM, ZCOM = self.COMdefine(self.x, self.y, self.z, np.s_[:]) #create a center of mass vector from all particles of type
        RCOM = (XCOM**2 + YCOM**2 + ZCOM**2) ** 0.5  #magnitude of com vector

        #determinging each particles position from first guess com
        x_ref = self.x - XCOM
        y_ref = self.y - YCOM
        z_ref = self.z - ZCOM

        RNEW = ( x_ref**2 + y_ref**2 + z_ref**2 ) ** 0.5 # absolute distance of each particle from COM
        RMAX = np.amax(RNEW) / 2 # set volume to half the maximum distance

        RCOM2 = RCOM - 1000 * u.kpc #initial change in com guess

        while np.absolute(RCOM - RCOM2) > delta * u.kpc: # recalculate com using half radii until convergence met
            RCOM = RCOM2
            vol_index = np.where(RNEW < RMAX) #finding which particles within new radii
            XCOM2, YCOM2, ZCOM2 = self.COMdefine( self.x, self.y, self.z, vol_index )
            RCOM2 = (XCOM2**2 + YCOM2**2 + ZCOM**2) ** 0.5
            RMAX = RMAX / 2 #reduce radius again
            x_ref2 = self.x - XCOM2
            y_ref2 = self.y - YCOM2
            z_ref2 = self.z - ZCOM2
            RNEW = ( x_ref**2 + y_ref**2 + z_ref2**2 ) ** 0.5
            #print(np.absolute(RCOM - RCOM2))

            # NOTE: 4) This iterative procedure is important, because when two galaxies collide
            # it will be difficult to distinguish which galaxy certain stars belong to.
            # Stars around the COM are likely to stay together the longest, so we need to
            # consider these stars especially when trying to determine the galaxies position.


        return XCOM2, YCOM2, ZCOM2

    def COM_V(self, XCOM, YCOM, ZCOM):
        # Determine COM velocity vector using only particles within 15 kpc of position COM

        x_ref = self.x - XCOM
        y_ref = self.y - YCOM
        z_ref = self.z - ZCOM

        RMAG = ( x_ref**2 + y_ref**2 + z_ref**2 ) ** 0.5 #absolute distance of each particle from COM

        vol_index = np.where(RMAG < 15.0 * u.kpc) #distance index

        VXCOM, VYCOM, VZCOM = self.COMdefine(self.vx, self.vy, self.vz, vol_index)

        return VXCOM, VYCOM, VZCOM






# EXAMPLE OF USING A CLASS
##########################
"""
# Create a Center of mass object for the MW
MWCOM = CenterOfMass("/home/agibbs/MW_000.txt", 2)

# Create a center of mass object for M31 and M33
M31COM = CenterOfMass("/home/agibbs/M31_000.txt", 2)
M33COM = CenterOfMass("/home/agibbs/M33_000.txt", 2)


# Calculate quantities for MW data
#MW_mass = MWCOM.total_mass()
#print("MW Disk Mass:", MW_mass)

# Calculate COM of MW, M31, and M33 with delta of 50 pc
MW_xcom, MW_ycom, MW_zcom = MWCOM.COM_P(0.1)
MW_vxcom, MW_vycom, MW_vzcom = MWCOM.COM_V(MW_xcom, MW_ycom, MW_zcom)

M31_xcom, M31_ycom, M31_zcom = M31COM.COM_P(0.1)
M31_vxcom, M31_vycom, M31_vzcom = M31COM.COM_V(M31_xcom, M31_ycom, M31_zcom)

M33_xcom, M33_ycom, M33_zcom = M33COM.COM_P(0.1)
M33_vxcom, M33_vycom, M33_vzcom = M33COM.COM_V(M33_xcom, M33_ycom, M33_zcom)

# Print COM for MW, M31, M33
print("\n MW COM:")
print('x:', MW_xcom, 'y:', MW_ycom, 'z:', MW_zcom)
print('vx:', MW_vxcom, 'vy:', MW_vycom, 'vz:', MW_vzcom)

print("\n M31 COM:")
print('x:', M31_xcom, 'y:', M31_ycom, 'z:', M31_zcom)
print('vx:', M31_vxcom, 'vy:', M31_vycom, 'vz:', M31_vzcom)

print("\n M33 COM:")
print('x:', M33_xcom, 'y:', M33_ycom, 'z:', M33_zcom)
print('vx:', M33_vxcom, 'vy:', M33_vycom, 'vz:', M33_vzcom)

# Calculate separations and relative velocities of galaxies
MW_M31_sep = ( (MW_xcom - M31_xcom)**2 + (MW_ycom - M31_ycom)**2 + (MW_zcom - M31_zcom)**2 )**0.5
MW_M31_vel = ( (MW_vxcom - M31_vxcom)**2 + (MW_vycom - M31_vycom)**2 + (MW_vzcom - M31_vzcom)**2 )**0.5

M31_M33_sep = ( (M31_xcom - M33_xcom)**2 + (M31_ycom - M33_ycom)**2 + (M31_zcom - M33_zcom)**2 )**0.5
M31_M33_vel = ( (M31_vxcom - M33_vxcom)**2 + (M31_vycom - M33_vycom)**2 + (M31_vzcom - M33_vzcom)**2 )**0.5

# Print separations and relative velocities
print('\n')

print("MW - M31 Separation:", MW_M31_sep)
print("MW - M31 Velocity:", MW_M31_vel)

print('\n')

print("M31 - M33 Separation:", M31_M33_sep)
print("M31 - M33 Velocity:", M31_M33_vel)

print('\n')
"""
