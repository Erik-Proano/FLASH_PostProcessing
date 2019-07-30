'''
  DESCRIPTION
    This script reads HDF5 output files from The FLASH code and
    extract its fields, data and computes the extent of the 
    mixing layer.

  AUTHOR
    Erik S. Proano
    Embry-Riddle Aeronautical University
'''

import yt
from yt.funcs import mylog
import h5py
import numpy as np
import scipy.interpolate
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot as plt
from matplotlib import cm, ticker
from matplotlib.ticker import MaxNLocator
from mpi4py import MPI
import os

#fname = 'cylindrical_rmi_2d_hdf5_chk_'
fname = "spherical_rmi_3d_hdf5_chk_"
location = os.getcwd()

# Initialize communicators
#comm = MPI.COMM.WORLD()
#Nproc = int(comm.Get_size())
#Pid = int(comm.Get_rank())
initfile = 0
finalfile = 150

r_max = (3.90, "cm")	# Max. Radius as a tuple with (Ma. Radius, "units")
p_res =  2*1024			# Number of desired radial bins 
r_tar = 2.5

# Divide equal amount of work to each worker (Done by master thread)
start = initfile
end = finalfile
nfiles = finalfile-initfile
#local_start = int(start + Pid*(finalfile-start)/Nproc)
#local_end = int(local_start + (finalfile-start)/Nproc)
#local_nfiles = local_end - local_start
#print('Processor ', Pid, 'will read from', local_start, 'to', local_end)
#comm.Barrier()
#t_start = MPI.Wtime()
time = np.array([])
rad_dens =np.zeros([p_res, nfiles])
#center = [0.0, 0.0, 0.5]
for i in range(start, end):
   if i > 9 and i <= 99:
     file = fname + "00" + str(i)
   elif i > 99 and i <= 999:
     file = fname + "0" + str(i)
   elif i > 999 :
     file = fname + str(i)
   else :
     file = fname + "000" + str(i)
   
   mylog.setLevel(40)
   ds = yt.load(file)
   print("Reading", file)
   # Create a 4cm-radius sphere
   center = ds.domain_left_edge#-(0.45,0.45,0.)
   sp = ds.sphere(center, r_max)
   profile = yt.create_profile(sp, 'radius', ['pressure'], n_bins=p_res,
                               units = {'radius': 'cm', "pressure": "dyn/cm**2"},
                               logs = {'radius': False, "pressure": True})
   # Transform the profile from a dictionary to a numpy array
   profVal = list(profile.field_data.values())
   for k in profVal:
	   d = k
   rad_dens[:,i] = d	
   time = np.append(time, float(ds.current_time))
rad = np.array(profile.x)/r_tar
X, Y = np.meshgrid(time*1.E06, rad)
levels=MaxNLocator(nbins=512).tick_values(rad_dens.min(),
                                         rad_dens.max())
#rbf = scipy.interpolate.Rbf(time, rad, rad_dens,function='linear')
#zi = rbf(X,Y)
#plt.imshow(zi, vmin=z.min(), vmax=z.max(), origin='lower',
#           extent=[x.min(), x.max(), y.min(), y.max()])
#plt.scatter()
fig = plt.figure()
ax = fig.add_subplot(111)
cf = ax.contourf(X,Y,rad_dens, 
                 levels=levels,
                 locator=ticker.LogLocator(), 
                 cmap=cm.binary)
#ax  = fig.add_subplot(111, projection='3d')
#ax.plot_surface(X, Y, rad_dens, cmap=cm.binary,
#                               linewidth=0, antialiased=False)
ax.set_ylabel(r'$r^{*}$',fontsize=16)
ax.set_xlabel(r'Time ($\mu s$)',fontsize=16)
cbar = fig.colorbar(cf, ax=ax)
fig.tight_layout()
fig.savefig("rt.png")
