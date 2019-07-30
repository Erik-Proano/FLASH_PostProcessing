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
from sympy.physics.vector import ReferenceFrame
from sympy.physics.vector import curl

def _vort(field, data):
    #R = ReferenceFrame('R')
    #velfield = np.array(data["velocity_x"]*R.x+data["velocity_y"]*R.y+data["velocity_z"]*R.z)
    #vort = curl(velfield, R)
    dv_dx = np.diff(data["velocity_y"])/0.0001
    du_dy = np.diff(data["velocity_x"])/0.0001 #/data["dy"][:-1]
    dv_dx = np.append(dv_dx, dv_dx[-1])
    #dv_dx = dv_dx/data["dx"]
    du_dy = np.append(du_dy, du_dy[-1])
    return dv_dx - du_dy

np = yt.numpy
fname = "cylindrical_rmi_2d_hdf5_chk_0005"
field = ("gas","vorticity_x")
res = 1024

ds = yt.load(fname)				# Load hdf5 output file
ds.add_field(("gas","vort"), function=_vort,
             sampling_type="cell", display_name="vorticity")
# Add num_ghost_zones=1 to add ghost cells for computing derivatives
cg = ds.covering_grid(level=0, left_edge=ds.domain_left_edge,
                      dims=ds.domain_dimensions,
                      num_ghost_zones=1)	# Create fixed covering grid
x = np.linspace(0., 4., res)
y = np.linspace(0., 4., res)
X, Y = np.meshgrid(x,y)
dens_grad_x, dens_grad_y = np.gradient(cg[field][:,:,0])
dens_grad_mag = np.sqrt(dens_grad_x**2+dens_grad_y**2)
slc = yt.plot_2d(ds, ("gas", "vort"))
#slc.set_log("vort", True)
slc.save()
#plt.contourf(X,Y,cg[field][:,:,0], 50)
#plt.set_cmap("binary")
#cbar = plt.colorbar()
#cbar.set_label(field[1]+" ["+str(cg[field].units)+"]")
#plt.show()
