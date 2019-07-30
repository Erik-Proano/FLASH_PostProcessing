# Parallel Post-Processing Interface for HDF5 Files from the Flash Code
Erik Proano
---
## Post Module
### getFileName
Returns an array containing all the file names based on the case name and the initial and final file indexes which are required as input arguments by this function
### ReynoldsDecomp
Performs a Reynolds decomposition in order to compute average flow and fluctuating quantities of all the data files
## Plot Functions
### plot_gradient
Saves in disk a surface plot with the gradient of a desired thermodynamic or kynematic flow quantity
### plot_ICs
Saves on disk surface plots of the initial conditions regarding thermodynamic and kinematic quantities of the flow i.e. the initialization of the flow within the Flash code framework
### plot_integralMixLayWidth
Returns data and saves plot in disk of the integral mixing layer width development in time of a cylindircal or spherical density interface
### plot_mixFraction
Computes statistical molecular fraction for quantifying mixing between two different species
### plot_MIxLayDevelop
Returns data and saves plot in disk of the real mixing layer width development in time of a cylindircal or spherical density interface. The data computed using this function contains noise compared to the data computed with the integral mixing layer width
### plot_shockConditions
Saves in disk plots of the conditions (ratios) ahead and behind an imploding (or exploding) cylindircal or spherical shock wave
### plot_vorticity
Saves in disk a surface plot of the vorticity field of a given simulation
### plot_zoom
Saves in disk zoomed surface plots of a desired variable
### read_hdf5_genSurfPlot
Saves in disk surface plots of the whole simulation hdf5 files given a desired flow variable
### read_hdf5_parallel
Computes the real mixing layer width using a desired number of processors. Execution time is greatly improved compared to the plot_MIxLayDevelop function. This function is also written in a more-general purpose manner; hence, allowing the user to easily edit the code depending on the geometry of the simulation and finally obtain faster results due to the implementation of MPI API for Python mpi4py and the yt API for reading hdf5 Flash output data

