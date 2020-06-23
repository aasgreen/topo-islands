import numpy as np
from numpy.random import rand
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import glob as glob
import trackpy as tp
import pandas as pd
import json
from scipy.interpolate import interp1d
from scipy import interp
from PIL import Image
from PIL import ImageFont
from PIL import ImageDraw
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from skimage import exposure,img_as_ubyte
from moviepy.editor import VideoClip
from moviepy.editor import ImageSequenceClip
from skimage import color
import datetime
import time
import argparse
import os
import h5py
from mpl_toolkits.mplot3d import Axes3D
plt.ion()
def read_files(args):
    '''read in, sort by time and output numpy arrays containing each image frame in the simulation'''
    
    #check if hdf5 file exists
    for fname in os.listdir(args.fName):
        if fname.endswith('h5'):
            if fname == 'defectGrid.h5':
                print('h5 defect file found')
                hf = h5py.File(args.fName+'/'+fname)
                fN = pd.DataFrame(list(hf.keys()))
                print(len(list(hf.keys())))
#            dfN = pd.DataFrame(glob.glob(args.fName+'/defect*.dat'))
#                params = np.loadtxt('params.txt')
#            print(fN)
                fN['time'] = fN[0].str.extract(r'(\d*\.?\d*)\.dat').astype('float')
                fN.sort_values(by=['time'],inplace=True)
#            dfN['time'] = dfN[0].str.extract(r'(\d*\.?\d*)\.dat').astype('float')
#            dfN.sort_values(by=['time'],inplace=True)
            #Sort fileNames by number

                imSeq = [np.array(hf.get(f)) for f in fN[0]]
                print(imSeq)
            break
 
    return [imSeq, fN]

if __name__==   "__main__":
    parser = argparse.ArgumentParser(description='Make some movies')
    parser.add_argument('fName',type=str,help='directory name for data')
    args = parser.parse_args()

    imSeq, fN = read_files(args)

   #first, find all positive defects using pandas
    test = imSeq[100] 
    posD = np.where(test >0)
    feature = { 'x': posD[0], 'y':posD[1], 'frame': 100}
    pfeatures = pd.DataFrame()
    nfeatures = pd.DataFrame()
    for i,frame in enumerate(imSeq):
        posD = np.where(frame > 0)
        negD = np.where(frame <0)
        #print(posD)
        if posD[0].size != 0:
            for x,y in zip(posD[0],posD[1]):
                pfeatures =pfeatures.append([{ 'x': x, 'y':y, 'frame': i}])
 #print(posD)
        if negD[0].size != 0:
            for x,y in zip(negD[0],negD[1]):
                nfeatures =nfeatures.append([{ 'x': x, 'y':y, 'frame': i}])


#now link these
pt = tp.link_df(pfeatures, 1,memory=4)
nt = tp.link_df(nfeatures, 1,memory=4)
#tp.plot_traj(pt)
#tp.plot_traj(nt)

#create worldline plot. This will plot the x, y positions of the defects, with the frame number being the z coordinate.

pz = pt['frame']
px = pt['x']
py = pt['y']


nz = nt['frame']
nx = nt['x']
ny = nt['y']

fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')
ax.scatter(px,py,pz, c='r', marker = '+', alpha = .7, s = 3)
ax.scatter(nx,ny,nz, c='b', marker = 'o', alpha = .7, s = 3)
