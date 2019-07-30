##
## @brief Plot CHisnell and simulation data shock conditions
##        This script is intended to validate the implementation of
##        Chisnell's approximation for cylindrical and spherical 
##        shock waves inside the FLASH code for the cases developed
## 
## @author Erik Proano
##

import os
import numpy as np
import matplotlib.pyplot as plt

## Set Chisnell data file and simulation shock data file
fileRoot = "/home2/proanoe/Documents/codes/shock_chisnell/fort.7"
fileSim  = os.getcwd() + "/fort.7"

## Load data file to numpy arrays
dataRoot = np.loadtxt(fileRoot)
dataSim  = np.loadtxt(fileSim)

## Organize Chisnell's data big array into individual 1d arrays
R_rat_C = dataRoot[:,0]
rho_rat_C = dataRoot[:,1]
u_rat_C = dataRoot[:,2]
pres_rat_C = dataRoot[:,3]

## Organize simulation's data array into individual 1d arrays
R_rat_S = dataSim[:,0]
rho_rat_S = dataSim[:,1]
u_rat_S = dataSim[:,2]
pres_rat_S = dataSim[:,3]

#rho_error = (np.abs(rho_rat_S-rho_rat_C)/rho_rat_C)*100
#u_error = (np.abs(u_rat_S-u_rat_C)/u_rat_C)*100
#pres_error = (np.abs(pres_rat_S-pres_rat_C)/pres_rat_C)*100

R_rat_S = 1/R_rat_S		# Radial data in comparable form as Chisnell

opt = 1
print("Plot Options")
print("1 - Density ratio")
print("2 - Velocity ratio")
print("3 - Pressure ratio")
print("0 - Exit script")
while (opt != 0):
	opt = input("Enter option: ")
	opt = int(opt)
	if (opt == 0):
		print("Exiting now....")
	elif (opt == 1):
		plt.plot(R_rat_C, 1/rho_rat_C,'b-')
		plt.plot(R_rat_S, 1/rho_rat_S,'r.')
		plt.xlim(R_rat_C.max(), R_rat_C.min())
		plt.xlabel(r"$\xi^{-1}$", fontsize = 14)
		yl = plt.ylabel(r"$\frac{\rho_s}{\rho}$", fontsize=14, labelpad=10)
		yl.set_rotation(0)
		plt.xlim([1.0, 0.5])
		plt.legend(["Chisnell data", "Simulation data"])
		plt.show()
	elif (opt == 2):
		plt.plot(R_rat_C, u_rat_C,'b-')
		plt.plot(R_rat_S, u_rat_S,'r.')
		plt.xlim(R_rat_C.max(), R_rat_C.min())
		plt.xlabel(r"$\xi^{-1}$", fontsize = 14)
		yl = plt.ylabel(r"$\frac{u_s}{u}$", fontsize=14, labelpad=10)
		yl.set_rotation(0)
		plt.xlim([1.0, 0.5])
		plt.legend(["Chisnell data", "Simulation data"])
		plt.show()
	elif (opt == 3):
		plt.plot(R_rat_C, pres_rat_C,'b-')
		plt.plot(R_rat_S, pres_rat_S,'r.')
		plt.xlim(R_rat_C.max(), R_rat_C.min())
		plt.xlabel(r"$\xi^{-1}$", fontsize = 14)
		yl = plt.ylabel(r"$\frac{p_s}{p}$", fontsize=14, labelpad=10)
		yl.set_rotation(0)
		plt.xlim([1.0, 0.5])
		plt.legend(["Chisnell data", "Simulation data"])
		plt.show()
	else:
		print("Invalid option, try again!")
	
