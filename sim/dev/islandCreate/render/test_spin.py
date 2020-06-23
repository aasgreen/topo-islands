import numpy as np
from numpy.random import rand
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import json
from skimage import exposure,img_as_ubyte
import time
import argparse
import os
import glob
import re
import h5py

def binMatrix(arr, new_shape):
    shape = (new_shape[0], arr.shape[0] // new_shape[0],
            new_shape[1], arr.shape[1] // new_shape[1])
    return arr.reshape(shape).mean(-1).mean(1)



parser = argparse.ArgumentParser(description='Make some movies')
parser.add_argument('fName',type=str,help='directory name for data')
args = parser.parse_args()

anchor = np.loadtxt(args.fName+'/anchor.dat')
hf= h5py.File(args.fName+'/anchor.h5')
anchor_h = np.array(hf.get('anchor.dat'))
hf.close()


spin = anchor%(2*np.pi)

spin = spin.T #read from fortran ascii, needs to be transposed
#check that spin and anchor_h are the same

#now, test to make sure the quiver plots are working correctly
Y,X = np.indices(spin.shape)
U = np.cos(spin)
V = np.sin(spin)

[X,Y,U,V] = (binMatrix(t,[50,50]) for t in [X,Y,U,V])

fig, ax = plt.subplots(figsize=(18,16))
ax.imshow(spin, cmap='twilight', alpha=.4,origin='top')
ax.quiver(X,Y,U,V, angles ='uv')

#the arrows should line up with the correct angle



