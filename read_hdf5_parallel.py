##  @package read_hdf5_parallel.py
#    @brief   This script reads HDF5 output files from The FLASH code and
#    	    extract its fields, data and computes the extent of the 
#           mixing layer by rotating a line buffer fixed at one end. The
#           buffer catches mass fraction along it and determines the 
#           extend of the mixing layer by a maximum minus a minimum value
#           in the turbulent mixing region denoted by mfrac!=1 0r mfrac !=0.
#
#    @author Erik S. Proano
#     	      Embry-Riddle Aeronautical University

import matplotlib
#matplotlib.use('agg')	 # If matplotlib error arises, uncomment this line
import matplotlib.pyplot as plt
import yt
from yt.funcs import mylog
import h5py
import numpy as np
from mpi4py import MPI
import os
import sys

# Establish base name and path for the simulation output files
fname = "cylindrical_rmi_2d_hdf5_chk_"
fieldName = 'mach_number'
location = os.getcwd()			# Get current directory path
#os.mkdir(fieldName + '_map')
min_field = 1.0E-03
max_field = 5E-02
## The code will read from file 'fname_initfile' to 'fname_finalfile'
initfile = 5				# Initial file index
finalfile = 100				# Final file index
nfiles = finalfile - initfile
p_MixLayr = np.array([])
p_time = np.array([])
yt.enable_parallelism()
# Initialize MPI Communicator
comm = MPI.COMM_WORLD
Nproc = int(comm.Get_size())
Pid = int(comm.Get_rank())
master = 0				# Define Rank 0 as master thread
# Divide equal amount of work to each worker (Done by master thread)
start = initfile
local_start = int(start + Pid*(finalfile-start)/Nproc)
local_end = int(local_start + (finalfile-start)/Nproc)
local_nfiles = local_end - local_start

## Display the correct message if multiple files are read/rank or just one/rank
if local_start != local_end - 1:
  print('Processor ', Pid, 'will read from index file', local_start, 'to', local_end-1)
else:
  print("Processor ", Pid, "will only read index file", local_start)
#comm.Barrier()
if Pid == master:
   print('Current Directory: ',location)
   print('Number of Processors: ', Nproc)
   print('\n')
   print('----------------File-----------------', \
         '---------------time(ms)--------------', \
         '---Processor---\n')
comm.Barrier()
t_start = MPI.Wtime()		# Start run time measurement
for i in range(local_start, local_end):    # complete the file name with the correct index
   if i > 9 and i <= 99:
     file = fname + "00" + str(i)
   elif i > 99 and i <= 999:
     file = fname + "0" + str(i)
   elif i > 999 :
     file = fname + str(i)
   else :
     file = fname + "000" + str(i)

   # Load file
   mylog.setLevel(40)                  # Show INFO only in command prompt
   ds = yt.load(file)
   #ad = ds.all_data()
   #fd = ad['gas', fieldName]
   print('   ',file,'\t\t', float(ds.current_time)*1000, \
         '\t\t', Pid)

   rad = 15.0
   res = 2048
   xmix = np.array([])
   ymix = np.array([])
   rmix = np.array([])
   mylog.setLevel(40)           	# Hide output log in command line
   for j in range(0, 90):		# Loop for changing the angle of the line from 0 to 90 deg
      xp = rad*np.cos(j*np.pi/180)		# Compute x coordinate of the end point
      yp = rad*np.sin(j*np.pi/180)	 	# Compute y coordinate of the end point
      start = [0.0, 0.0, 0.5]
      end = [xp, yp, 0.5]
      lb = yt.LineBuffer(ds, start, end, res)	# Generate the line and get data acroos line
      #srt = np.argsort(ray['x'])
      xc = np.array(lb['gas', 'x'])
      yc = np.array(lb['gas', 'y'])
      mfrac = np.array(lb['fld1'])
      for k in range(0, len(mfrac)):	            # Start searching in the line (ray)
         if mfrac[k] >= 0.01 and mfrac[k] <= 0.98:  # Save coordinates of the mixture
           xmix = np.append(xmix, xc[k])
           ymix = np.append(ymix, yc[k])
           rmix = np.append(rmix, np.sqrt(xmix[len(xmix)-1]**2 + ymix[len(ymix)-1]**2))
           #p = p + 1
             
   p_MixLayr = np.append(p_MixLayr, max(rmix) - min(rmix))
   p_time = np.append(p_time, float(ds.current_time))
   #plt.plot(xmix, ymix, '.')
   #plt.savefig(location + '/' + fieldName + '_map/' + file + '_MixingLayer.png')
   #def plot_fields(fieldName):              # Generate and save density maps
     #p = yt.plot_2d(ds, 'z', fieldName)
     #p.set_log(fieldName, False)
     #p.set_cmap(fieldName, 'RdBu_r')
     #slc.set_zlim('dens', min_field, max_field)
     #slc.set_axes_unit('cm')
     #slc.hide_axes()
     #slc.hide_colorbar()
     #p.save()
print("Processor %d computation terminated" % Pid)
comm.Barrier()
## Gather results from all ranks to an auxiliary variable
#  The aux variable size is (Nproc, nfiles)
MixLayr_aux = comm.gather(p_MixLayr, root=master)
time_aux = comm.gather(p_time, root=master)
if Pid == master:                # Everything is collected and organized by the Master
  print("All ranks are done with computations!")
  MixLayr = np.array([])	# Initialize mix. layer array
  time = np.array([])		# Initialize time array
  ## Append values from aux. variable to final 1D array
  for i in range(0, np.array(MixLayr_aux).shape[0]):
     MixLayr = np.append(MixLayr, MixLayr_aux[i])
     time = np.append(time, time_aux[i])
  #MixLayr = np.array([])
  #time = np.array([])
  #comm.Gather(p_MixLayr, MixLayr, root=0)
  #comm.Gather(p_time, time, root=0)     
  #for i in range(0, Nproc):            #thread
    #comm.recv(p_MixLayr, source=i, tag=77)
    #comm.recv(p_time, source=i, tag=77)
    #print(np.array([p_time]).T, np.array([p_MixLayr]).T) 
    #for j in range(0, len(p_MixLayr)):
      #MixLayr = np.append(MixLayr, p_MixLayr[j])
      #time = np.append(time, p_time[j])
     
t_end = MPI.Wtime()		# End run time measurement

if Pid == 0:		## Master rank generates all the results
  runtime = t_end - t_start
  print('Total runtime is', runtime, 'sec')
  print(np.array([time]).T, np.array([MixLayr]).T) 
  f = open('mixing_layer_width.dat', 'w')
  f.write('#  Time(s)\t\tWidth(cm)\n')
  for i in range(0, len(MixLayr)):
    f.write(str(time[i]) + '       ' + str(MixLayr[i]) + '\n')
  f.close()
  #plt.plot(1000*time, MixLayr, '.-')
  #plt.xlabel('time (ms)')
  #plt.ylabel('Mixing Layer Width (cm)')
  #plt.savefig('mixing_layer_width.png')
  len(time)
  print('Array size is ', np.array(time).shape)
  print("Mixing layer data placed in %s" % location)
  sys.exit("Execution terminated succesfully!") # Terminate with sys message
