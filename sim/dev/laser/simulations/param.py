import numpy as np
import random
#arguments go like g, beta, N, endT

numRuns = 1
kappaL = [i for i in np.linspace(1,1,1)]
betaL = [i for i in np.linspace(2,2,1)]
seed = [int(random.randint(0,32000)) for i in np.arange(numRuns)]
kappa = 1
mu = 0.
beta = 1
betaLaser = .01
gridSize =100
endTime = 200
s = seed[0]
#params = np.array([(k,beta,mu, gridSize,endTime,s, betaLaser) for k in kappaL for beta in betaL for s in seed])
params = np.array([[kappa, beta, mu, gridSize, endTime, s, betaLaser]])
print(params)
print(kappaL)
np.savetxt('params.txt',params,delimiter=' ',fmt=['%03.2f','%04.3f','%04.3f','%d','%d','%d', '%03.2f'])


