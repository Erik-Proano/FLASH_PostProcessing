import postModules as post
import os
import numpy as np
from yt.funcs import mylog
import yt
from matplotlib import pyplot as plt

fname = "cylindrical_rmi_2d_hdf5_chk_"	   # Establish base name and path for the simulation output files
#fieldName = 'mach_number'
fieldName = ['density']
dirName = fieldName[0] + '_map'
location = os.getcwd()
fullDirPath = location + '/' + dirName
path = "/home2/proanoe/Simulations2018/2D/single_mode/run1_Ding_case1/quarter/a_1mm/test6/"

# mproc.cpu_count() <= Nproc
# Initialize MPI Communicator
#comm = MPI.COMM_WORLD
#Nproc = int(comm.Get_size())
#Pid = int(comm.Get_rank())
Pid=0
Nproc=1

    
min_field = 1.0E-03
max_field = 5E-02
initfile = 0
finalfile = 50
nfiles = finalfile - initfile
nSnaps = nfiles
nVars = len(fieldName)		# number of variables composing the variable
res = 1024					# Get resolution of a squared-cartesian grid
#q = np.zeros([res,res,nvars,nSnaps])	# Initialize variable for correlation
b = np.zeros([nfiles])
time = np.zeros([nfiles])

# Divide equal amount of work to each worker (Done by master thread)
start = initfile
local_start = int(start + Pid*(finalfile-start)/Nproc)
local_end = int(local_start + (finalfile-start)/Nproc)
local_nfiles = local_end - local_start
if Pid == 0:
	files = post.getFileName(nfiles, path, fname)
	average = post.ReynoldsDecomp(nfiles, path, fieldName, files)

for i in range(local_start, local_end):    # complete the file name with the correct index
   print( "Reading file ", files[i], ' with processor ', Pid)
   # Load file
   mylog.setLevel(40)                  # Show no INFO in command prompt
   ds = yt.load(path+files[i])
   #all_data = ds.covering_grid(0, ds.domain_left_edge, [1024,1024,1])	# Extract data in readable form
   #all_data = ds.all_data()
   dens = np.array(ds.r['density'])
   dens_f = dens-average[:,0]
   b[i] = np.average(dens_f*dens_f)/average[:,0]
   time[i] = float(ds.current_time)
   
plt.plot(time*1.E-03, b)
plt.xlabel('Time (ms)')
plt.ylabel('b')
plt.grid()
plt.savefig("./bCorrelation.png")
   
