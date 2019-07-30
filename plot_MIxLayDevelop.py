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

import numpy as np
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter

def read_datFile(fname, time, width, time_gr, growth):
  fdata = np.loadtxt(fname)
  time = fdata[:,0]
  width = fdata[:,1]
  size = len(width)
  #growth = np.zeros(size-1)
  #time_gr = np.zeros(size-1)
  accel = np.zeros(size-2)
  time_ac = np.zeros(size-2)
  for i in range(1,size):
     growth[i-1] = (width[i]-width[i-1])/(time[i]-time[i-1])
     time_gr[i-1] = time[i]
     if i <= size-2:
        accel[i-1] = (width[i]-2*width[i+1]+width[i-1])/(time[i]-time[i-1])
  return time, width, time_gr, growth

fname = 'mixing_layer_width.dat'
#location = ['/home2/proanoe/Simulations2018/2d_multimode/wavenumber_effects/a1_k16_6/',\
#            '/home2/proanoe/Simulations2018/2d_multimode/wavenumber_effects/a1_k6_32/',\
#            '/home2/proanoe/Simulations2018/2d_multimode/wavenumber_effects/a1_k6_64/',\
#            '/home2/proanoe/Simulations2018/2d_multimode/wavenumber_effects/a1_k16_64/',\
#            "/home2/proanoe/Simulations2018/2d_multimode/wavenumber_effects/a1_k32_64/",\
#            "/home2/proanoe/Simulations2018/2d_multimode/wavenumber_effects/a1_k6_16_32_64/"]
#legend = [r'$k_1=6, k_2=16$','$k_1=6, k_2=32$','$k_1=6, k_2=64$', 
#           "$k_1=16, k_2=64$", "$k_1=32, k_2=64$", "$k=all$"]
location = ["/home2/proanoe/Simulations2018/2D/single_mode/a1_k6/", \
            "/home2/proanoe/Simulations2018/2D/single_mode/a1_k16/", \
            "/home2/proanoe/Simulations2018/2D/single_mode/a1_k32/", \
            "/home2/proanoe/Simulations2018/2D/single_mode/a1_k64/a1_k64/"]
legend = ["k=6","k=16", "k=32","k=64"]
log_flag = False
fdata = np.loadtxt(location[0]+fname)
time = fdata[:,0]
width = fdata[:,1]
size = len(width)
growth = np.zeros(size-1)
time_gr = np.zeros(size-1)
lc = ['b','r','k',"g","y", "m"]
for i in range(0,len(location)):
  read_datFile(location[i]+fname, time, width, time_gr, growth)
  # Filter data for smoothing
  tt = np.linspace(time_gr.min(), time_gr.max(), 1000)
  itp = (interp1d(time_gr, growth, kind='linear'))
  growth_fil = savgol_filter(itp(tt), 101, 3)
  # Plot data
  #plt.plot(time_gr*1E06,growth/100)
  plt.figure(1)
  plt.plot(tt*1E6, growth_fil/100, lc[i])
y = 0.01*((tt[:20])**-0.5)
#plt.plot(time[:20]*1E06, y, "k--")
#plt.annotate(r"~$\mathbf{t^{-0.5}}$", xy = (75,2.5))
plt.xlabel(r'time ($\mu s$)')
plt.ylabel(r'$\dot{a}$ (m/s)')
plt.grid(True)
axes = plt.gca()
if log_flag:
	axes.set_xscale("log")
	axes.set_yscale("log")
#axes.set_xlim([0, 300])
plt.legend(legend, loc='lower left')
plt.show()
for i in range(0,len(location)):
  #read_datFile(location[i]+fname, time, width, time_gr, growth)
  fdata = np.loadtxt(location[i]+fname)
  time = fdata[:,0]
  width = fdata[:,1]
  plt.figure(2)
  plt.plot(time*1E06, width/2.5, lc[i])
y = 30*((time[:20])**0.5)
#plt.plot(time[:20]*1E06, y/2.5, "k--")
#plt.annotate(r"~$\mathbf{t^{0.5}}$", xy = (125,0.3))
plt.xlabel(r'time ($\mu s$)')
plt.ylabel(r'$a/R_{0}$')
plt.grid(True)
axes = plt.gca()
if log_flag:
	axes.set_xscale("log")
	axes.set_yscale("log")
#axes.set_xlim([0, 300])
plt.legend(legend, loc='upper left')
plt.show()

