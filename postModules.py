## @package postModules
#  @brief This module contain useful functions for post-processing
#         simulation data from the FLASH code.
#
#  Functions provided in the present module describe general
#  processes used commonly for the analysis of turbulent mixing
#  @author Erik Proano
import numpy as np
from yt.funcs import mylog
from pathlib import Path
import yt

## getFileName return a tuple with the file names
#
#  Based on a base name, this method returns the complene name of files
#  depending on the number of files requested.
#  @param nfiles    The total snapshot count
#  @param path      Path to directory cantaining files
#  @param basename  Common name for the files
def getFileName(initfile, finalfile, path, basename):
	nfiles = finalfile - initfile
	files = [None]*nfiles
	for i in range(initfile, finalfile):
		# Append the correct index to each file
		if i > 9 and i <= 99:
			file = basename + "00" + str(i)
		elif i > 99 and i <= 999:
			file = basename + "0" + str(i)
		elif i > 999 :
		    file = basename + str(i)
		else :
			file = basename + "000" + str(i)
		# Check if file exist or path is correct		
		try:
			my_abs_path = Path(path+file).resolve()
		except FileNotFoundError:
			print("{} file does not exist or path is wrong!\n".format(file))
		else:
			mylog.setLevel(40)
			# print("Reading file {}".format(file))
			files[i-initfile] = file
	return files

## This function performs a Reynolds Decomposition
#  
#  The function extracts the mean flow and the fluc-
#  tuating part of the flow for turbulent flow simu-
#  lations carried out by the FLASH code.
#  @param nfiles    The total snapshot count
#  @param path      Path to directory containing files
#  @param fieldname List with the data fields of interest
#  @param basename  Shared name for the files
def ReynoldsDecomp(nfiles, path, fieldname, files):
	print("Starting with averaging process...")
	time = np.zeros(nfiles)				# Initialize time array
	#files = getFileName(nfiles, path, basename)
	j = 0
	for k in range(0,len(fieldname)):
		print("Starting averaging of "+ fieldname[k])
		for i in files:
			print(i)		   
			ds = yt.load(path+i)									# Load the file to the scope
			data = ds.r[fieldname[k]]
			time[j] = ds.current_time
			if j == 0:
				averaged_data = np.zeros(len(data))					# Initialize accumulator
				if k == 0:
					average = np.zeros((len(data),len(fieldname)))	# Initialize array
			j += 1
			averaged_data = averaged_data + data
		average[:,k] = averaged_data/nfiles
		j = 0
	print("Done with averaging process...")
	return average
