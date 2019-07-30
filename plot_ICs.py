###################################################################
##  DESCRIPTION
##    This script reads HDF5 output files from The FLASH code and
##    extract its fields, data and computes the extent of the 
##    mixing layer.
##
##  AUTHOR
##    Erik S. Proano
##    Embry-Riddle Aeronautical University
##
###################################################################

import yt
from yt.funcs import mylog
import h5py
import numpy as np
from matplotlib import pyplot as plt
from mpi4py import MPI
import os

file = 'cylindrical_rmi_2d_hdf5_chk_0000'
figFormat = 'png'
location = os.getcwd()
print(location)
field = 'velocity_magnitude'
rad_0 = 2.5
var_0 = 34029
ds = yt.load(file)
rad_max = 4
xp = rad_max
yp = rad_max
start = [0.0, 0.0, 0.5]
end = [xp, yp, 0.5]
lb = yt.LineBuffer(ds, start, end, 2048)
xc = np.array(lb['gas', 'x'])
yc = np.array(lb['gas', 'y'])
var_init = np.array(lb['gas', field])
rad = np.sqrt(xc**2 + yc**2)
plt.plot(rad/rad_0, -var_init/var_0)
plt.xlabel(r'$\frac{R}{R_0}$')
plt.ylabel(r'$\frac{u_r}{a_0}$')
plt.grid()
plt.savefig(location + '/' + 'ICs_' + field + '.' + figFormat)
