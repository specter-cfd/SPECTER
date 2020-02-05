import numpy as np
import matplotlib.pyplot as plt

# Reads a binary file and plots a cut in the x-z plane.
# Execute in a shell with 'python plot_binary.py'
# Execute in ipython with '%run plot_binary.py'


# Path to the run data
path = '../../bin/'

# Grid data
x, y, z = [ np.loadtxt(path+i+'.txt') for i in ['x','y','z'] ]
nx, ny, nz = x.size, y.size, z.size
shape = (nx, ny, nz)

# Reads binary files
# For simulations in single precision change np.float32
# to np.float64
vx = np.fromfile(path+'out/vx.0001.out',dtype=np.float32).reshape(shape,order='F')

# Show a vertical cut of the field in the middle of the box
plt.figure()
plt.imshow(vx[:,ny//2,:].T, origin='lower', extents=[ x[0],x[-1],z[0],z[-1] ],
           interpolation='gaussian')
plt.xlabel('x')
plt.ylabel('z')
plt.colorbar()
plt.show()
