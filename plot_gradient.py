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
fname = "cylindrical_rmi_2d_hdf5_chk_0005"
field = ("gas","pressure")
res = 1024

ds = yt.load(fname)				# Load hdf5 output file
# Add num_ghost_zones=1 to add ghost cells for computing derivatives
cg = ds.covering_grid(level=0, left_edge=ds.domain_left_edge,
                      dims=ds.domain_dimensions,
                      num_ghost_zones=1)	# Create fixed covering grid
x = np.linspace(0., 4., res)
y = np.linspace(0., 4., res)
X, Y = np.meshgrid(x,y)
field_grad_x, field_grad_y = np.gradient(cg[field][:,:,0])
field_grad_mag = np.sqrt(field_grad_x**2+field_grad_y**2)
plt.contourf(X,Y,field_grad_mag, 50)
plt.set_cmap("binary")
cbar = plt.colorbar()
cbar.set_label(field[1]+" ["+str(cg[field].units)+"]")
plt.show()
