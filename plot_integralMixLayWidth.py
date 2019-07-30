#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  plot_integralMixLayWidth.py
#  
#  Copyright 2019 Erik Proano <erik.proano@protonmail.com>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  
#  @brief Compute and plot the molecular mixture fraction
#
#  @author  Erik Proano
#


import yt
import os
from yt.funcs import mylog
import postModules as post
from mpi4py import MPI
import matplotlib.pyplot as plt
from yt import numpy as np

def _mixture(field, data):
    return data["fld1"]*(1-data["fld1"])

#fname = "cylindrical_rmi_2d_hdf5_chk_"	# 2D
fname = "spherical_rmi_3d_hdf5_chk_"	# 3D
location = os.getcwd()
path = location + '/'

yt.add_field(("gas","mix"), function=_mixture, take_log=False,
              units="", display_name="Mix Fraction Numerator", 
              sampling_type="cell")

r_max = (3.90, "cm")	# Max. Radius as a tuple with (Ma. Radius, "units")
p_res =  1*1024			# Number of desired radial bins 
r_tar = 2.5

# Initialize MPI Communicator
comm = MPI.COMM_WORLD
Nproc = int(comm.Get_size())
Pid = int(comm.Get_rank())
master = 0
if Pid == 0:
	initfile = int(input("Enter first output file index: "))
	finalfile = int(input("Enter last output file index: "))
	nfiles = finalfile - initfile
else:
	finalfile = 0
	initfile  = 0
	nfiles = 0
initfile = comm.bcast(initfile, root=0)
finalfile = comm.bcast(finalfile, root=0)
nfiles = comm.bcast(nfiles, root=0)
mfrac1 = np.zeros([p_res, nfiles])

# Divide equal amount of work to each worker (Done by master thread)
start = initfile
local_start = int(start + Pid*(finalfile-start)/Nproc)
local_end = int(local_start + (finalfile-start)/Nproc)
local_nfiles = local_end - local_start

if Pid == master:
    files = post.getFileName(initfile, finalfile, path, fname)
    time = np.zeros(nfiles)
    mixFrac = np.zeros(nfiles)
else:
    #files, time, mixFrac = None, None, None
    files = None
    time = None
    mixFrac = None
p_mixFrac = np.zeros(local_nfiles)
p_time = np.zeros(local_nfiles)
comm.Barrier()
## Display the correct message if multiple files are read/rank or just one/rank
if local_start != local_end - 1:
  print('Processor ', Pid, 'will read from index file', local_start, 'to', local_end-1)
else:
  print("Processor ", Pid, "will only read index file", local_start)

files = comm.bcast(files, root=master)
for i in range(local_start, local_end):
   print( "Reading file ", files[i], ' with processor ', Pid)
   mylog.setLevel(40)
   ds = yt.load(path+files[i])
   # Create a 4cm-radius sphere
   center = ds.domain_left_edge
   sp = ds.sphere(center, r_max)
   profile = yt.create_profile(sp, 'radius', ['fld1','mix'], n_bins=p_res,
                               units = {'radius': 'cm', "fld1": "", "mix": ""},
                               logs = {'radius': False, "fld1": False, "mix": False})
   mfrac1[:,i] = np.array(profile["fld1"])	
   p_time[i-local_start] = float(ds.current_time)
   rad = np.array(profile.x)
   num = np.array(profile["mix"])
   num = np.trapz(num, rad)
   den = 4*mfrac1[:,i]*(1-mfrac1[:,i])
   den = np.trapz(den, rad)
   p_mixFrac[i-local_start] = den
comm.barrier()						# Wait for all
comm.Gather(p_time, time, root=master)	# Gather time in one array
comm.Gather(p_mixFrac, mixFrac, root=master)			# Gather enstrophy in one array
MPI.Finalize()						# Finalize communication

if Pid == master:
    f = open('mixing_zone_width_integral.dat', 'w')
    f.write('#  Time(s)\t\tMixing Fraction\n')
    for i in range(0, len(mixFrac)):
        f.write(str(time[i]) + '       ' + str(mixFrac[i]) + '\n')
    f.close()

    plt.plot(time*1.E+06, mixFrac/r_tar, "b")
    plt.xlabel(r"Time $\left(\mu s\right)$")
    #plt.ylabel(r"$\delta$")
    plt.ylabel(r"$a_0/R_0$")
    plt.savefig("mixing_zone_width_integral.png")

