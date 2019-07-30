####################################################
##
##  @brief  Plot gradient of a field
##      This script provides the gradient field
##     of a FLASH native field and plot it.
##
##     @author  Erik Proano
##
#################################################### 

import yt
import matplotlib.pyplot as plt

np = yt.numpy
fname = "cylindrical_rmi_2d_hdf5_chk_0050"
field1 = ("gas","pressure")
field2 = ("gas","density")

ds = yt.load(fname)				# Load hdf5 output file
res = ds.domain_dimensions[0]
# Add num_ghost_zones=1 to add ghost cells for computing derivatives
cg = ds.covering_grid(level=0, left_edge=ds.domain_left_edge,
                      dims=ds.domain_dimensions,
                      num_ghost_zones=1)	# Create fixed covering grid
x = np.linspace(0., 4., res)
y = np.linspace(0., 4., res)
X, Y = np.meshgrid(x,y)
pres_grad_x, pres_grad_y = np.gradient(cg[field1][:,:,0])
dens_grad_x, dens_grad_y = np.gradient(cg[field2][:,:,0])
pres_grad_mag = np.sqrt(pres_grad_x**2+pres_grad_y**2)
dens_grad_mag = np.sqrt(dens_grad_x**2+dens_grad_y**2)
vortMis = np.matmul(pres_grad_mag, dens_grad_mag)
plt.contourf(X,Y,vortMis, 50)
plt.set_cmap("binary")
cbar = plt.colorbar()
cbar.set_label("Product"+" ["+str(cg[field1].units*cg[field2].units)+"]")
plt.show()
