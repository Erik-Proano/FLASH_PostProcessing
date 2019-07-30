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

#import matplotlib
#matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib import cm
import yt
from yt.funcs import mylog
import h5py
import numpy as np
#from matplotlib import pyplot as plt
from mpi4py import MPI
import os

#def _entropy(field, data):
#	return data["pressure"]/(data["density"]**data["gam"])

def plot_fields(ds, fieldName, path):     
  # Generate and save user-defined field surface plots
  p = yt.plot_2d(ds, fieldName)
  p.set_log(fieldName, False)
  p.set_cmap(fieldName, cm.binary)
  #slc.set_zlim('dens', min_field, max_field)
  #slc.set_axes_unit('cm')
  #slc.hide_axes()
  #slc.hide_colorbar()
  p.save(path)



#initfile = np.int(input("Enter first output file index: "))
#finalfile = np.int(input("Enter last output file index: "))
#nfiles = finalfile - initfile

# Initialize MPI Communicator
comm = MPI.COMM_WORLD
Nproc = int(comm.Get_size())
Pid = int(comm.Get_rank())

if Pid == 0:
  # Establish base name and path for the simulation output files
  fname = "cylindrical_rmi_2d_hdf5_chk_"
  #fieldName = 'mach_number'
  fieldName = input("Enter field name:")		# Input field
  dirName = fieldName + '_map'
  location = os.getcwd()
  fullDirPath = location + '/' + dirName
  print(location)
  print('Number of Processors requested: ', Nproc)
  dir_exist = os.path.isdir(fullDirPath)
  if dir_exist:
    print('Currently directory ' + fullDirPath + ' already exist!')
  else:
    print('Creating directory ' + dirName + ' in' + location)
    os.mkdir(dirName)
    print('Done creating directory ' + fullDirPath)
  initfile = int(input("Enter first output file index: "))
  finalfile = int(input("Enter last output file index: "))
  nfiles = finalfile-initfile
else:
  # Establish base name and path for the simulation output files
  fname = "cylindrical_rmi_2d_hdf5_chk_"
  #fieldName = 'mach_number'
  fieldName = ""
  dirName = ""#fieldName + '_map'
  location = os.getcwd()
  fullDirPath = ""#location + '/' + dirName
  initfile = 0
  finalfile = 0
  nfiles = 0
comm.Barrier()
fieldName = comm.bcast(fieldName, root=0)
dirName = comm.bcast(dirName, root=0)
location = comm.bcast(location, root=0)
fullDirPath = comm.bcast(fullDirPath, root=0)
initfile = comm.bcast(initfile, root=0)
finalfile = comm.bcast(finalfile, root=0)
nfiles = comm.bcast(nfiles, root=0)
min_field = 1.0E-03
max_field = 5E-02

# Divide equal amount of work to each worker (Done by master thread)
start = initfile
local_start = int(start + Pid*(finalfile-start)/Nproc)
local_end = int(local_start + (finalfile-start)/Nproc)
local_nfiles = local_end - local_start
for i in range(local_start, local_end):    # complete the file name with the correct index
   if i > 9 and i <= 99:
     file = fname + "00" + str(i)
   elif i > 99 and i <= 999:
     file = fname + "0" + str(i)
   elif i > 999 :
     file = fname + str(i)
   else :
     file = fname + "000" + str(i)

   print( "Reading file ", file, ' with processor ', Pid)
   # Load file
   mylog.setLevel(40)                  # Show no INFO in command prompt
   ds = yt.load(file)
   #ad = ds.all_data()
   #fd = ad['gas', fieldName]
            
   plot_fields(ds, fieldName, fullDirPath)
