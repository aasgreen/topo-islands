#Testing the fortran module

import langanLib
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
print(langanLib.langan.__doc__)

#parameters (will be moved to yml file)
n = 100
kappa = .1
mu = .01
beta = 2
lnoise = np.sqrt(12.0*1.0/beta)

#Testing Initialize Grid

grid = np.zeros([n,n], order='f')

langanLib.langan.init_random_seed()
langanLib.langan.initializegrid(grid,n)


#Testing Update Scheme

ngrid = np.zeros([n,n], order = 'f')
hgrid = np.zeros([n,n], order = 'f')
[ngrid0, hgrid0] = langanLib.langan.update(grid, lnoise, kappa, mu, n)

lgrid = []
ngrid = grid

#measure energy
hgrid2 = np.zeros([n,n], order = 'f')
for i in np.arange(n):
    for j in np.arange(n):
        hgrid2[i,j] = langanLib.langan.hamxy(i+1,j+1,ngrid0,kappa,mu,np.pi)

plt.imshow(hgrid2[20:80,20:80].T, origin = 'top')
plt.show()
#for t in np.linspace(0,100):
#    grid = ngrid
#    [ngrid, hgrid] =  langanLib.langan.update(grid, lnoise, kappa, mu, n)
#    lgrid.append(grid)
    
    
fig,ax = plt.subplots(nrows = 1, ncols = 3, figsize = (14,13))
ax[0].imshow(grid.T, cmap = 'twilight', origin = 'top')
ax[1].imshow(ngrid.T, cmap = 'twilight', origin = 'top')
ax[2].imshow(hgrid.T, origin = 'top')
plt.tight_layout()
plt.show()
