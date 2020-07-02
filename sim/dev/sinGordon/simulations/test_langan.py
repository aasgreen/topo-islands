#Testing the fortran module

import langanLib
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
print(langanLib.langan.__doc__)

#Testing Initialize Grid

n = 100
grid = np.zeros([n,n], order='f')

langanLib.langan.initializegrid(grid,n)

print(grid)
print(grid.shape)
plt.imshow(grid.T, cmap = 'twilight', origin = 'top')
plt.show()

#Testing Update Scheme

ngrid = np.zeros([n,n], order = 'f')
hgrid = np.zeros([n,n], order = 'f')

[ngrid, hgrid] = langanLib.langan.update(grid, .01, .1, .5, 100)

fig,ax = plt.subplots(nrows = 3, ncols = 1)
ax[0].imshow(grid.T, cmap = 'twilight', origin = 'top')
ax[1].imshow(ngrid.T, cmap = 'twilight', origin = 'top')
ax[2].imshow(hgrid.T, cmap = 'twilight', origin = 'top')
plt.show()
