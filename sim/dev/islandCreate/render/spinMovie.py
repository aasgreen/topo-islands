import numpy as np
from numpy.random import rand
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import glob as glob

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
from moviepy.video.io.bindings import mplfig_to_npimage

from skimage import color
import datetime
import time
import argparse
import os
import h5py
import re
def read_files(args):
    '''read in, sort by time and output numpy arrays containing each image frame in the simulation'''
    
    #check if hdf5 file exists
    for fname in os.listdir(args.fName):
        if fname.endswith('h5'):
            if fname == 'dsetf.h5':
                print('h5 file found')
                hf = h5py.File(args.fName+'/'+fname)
                fN = pd.DataFrame(list(hf.keys()))
                print(len(list(hf.keys())))
#            dfN = pd.DataFrame(glob.glob(args.fName+'/defect*.dat'))
                params = np.loadtxt('params.txt')
#            print(fN)
                fN['time'] = fN[0].str.extract(r'(\d*\.?\d*)\.dat').astype('float')
                fN.sort_values(by=['time'],inplace=True)
#            dfN['time'] = dfN[0].str.extract(r'(\d*\.?\d*)\.dat').astype('float')
#            dfN.sort_values(by=['time'],inplace=True)
            #Sort fileNames by number

                imSeq = [np.array(hf.get(f)) for f in fN[0]]
                break
#            dimSeq = [np.loadtxt(f) for f in dfN[0]]
    else:

        fileNames = glob.glob(args.fName+'/out*.dat')
#        dfileNames = glob.glob(args.fName+'/defect*.dat')
        fN = pd.DataFrame(fileNames)
        fN['time'] = fN[0].str.extract(r'(\d*\.?\d*)\.dat').astype('float')
        fN.sort_values(by=['time'],inplace=True)

#        dfN = pd.DataFrame(dfileNames)
#        dfN['time'] = dfN[0].str.extract(r'(\d*\.?\d*)\.dat').astype('float')
#        dfN.sort_values(by=['time'],inplace=True)
        #Sort fileNames by number

        imSeq = [np.loadtxt(f) for f in fN[0]]
#        dimSeq = [np.loadtxt(f) for f in dfN[0]]
#    return [imSeq,dimSeq]
    return [imSeq, fN]

def binMatrix(arr, new_shape):
    shape = (new_shape[0], arr.shape[0] // new_shape[0],
            new_shape[1], arr.shape[1] // new_shape[1])
    return arr.reshape(shape).mean(-1).mean(1)




parser = argparse.ArgumentParser(description='Make some movies')
parser.add_argument('fName',type=str,help='directory name for data')
args = parser.parse_args()

imSeq, fN = read_files(args)
tPattern = re.compile('.*beta-(\d*\.?\d*)')
beta = float(tPattern.search(args.fName)[1])
    
fps=24
secondsPframe_in = fN['time'].diff().mean()
fps_in = secondsPframe_in**-1

params = np.loadtxt('param.txt')
radius = np.loadtxt(args.fName+'/radius.dat')
anchor = np.loadtxt(args.fName+'./anchor.dat').T #from fortran ascii, need to transpose
f0 = imSeq[0]
imShape = f0.shape
Y,X = np.indices(imShape)
U = np.cos(f0)
V = np.sin(f0)

[X,Y,U,V] = (binMatrix(t,[50,50]) for t in [X,Y,U,V])

fig,[ax,ax2] = plt.subplots(figsize=(18,16),nrows=1,ncols=2)
ax.set_aspect(aspect='equal')
ax.axis('off')
plt.tight_layout()
im2=ax.imshow((anchor%(2*np.pi))*180./np.pi, cmap='twilight', alpha = .3, origin = 'top', vmin=0, vmax=360)
Q=ax.quiver(X,Y, U,V,scale=30)
im1=ax2.imshow( (imSeq[0]%(2*np.pi))*180./np.pi, cmap = 'twilight', origin = 'top', vmin =0, vmax = 360)
fig.colorbar(im1,ax=[ax,ax2], orientation = 'horizontal')

def schler(state):
    '''take given state and return schieren state
    '''
    return np.sin(2*state)**2
def make_frame(t):
    i = int(t*fps_in)
    img = schler(imSeq[i])

    imgSize = img.shape
    xx, yy = np.indices(imgSize)
    cx, cy = [i//2 for i in imgSize]
    r= np.sqrt((xx-cx)**2+(yy-cy)**2)
    return color.grey2rgb(np.asarray(im))

def make_dframe(t):
    i = int(t*fps)
    #img = schler(imSeq[i])
    #img = dimSeq[i]
    spins = binMatrix(imSeq[i],[50,50])
    imgSize = spins.shape
    #print(imgSize)

    U = np.cos(spins)
    V = np.sin(spins)

    Q.set_UVC(U,V)
    ax2.clear()
    im2=ax2.imshow( (imSeq[i]%(2*np.pi))*180./np.pi, cmap = 'twilight', origin = 'top', vmin = 0, vmax=360)
    fig.suptitle('t: {:04.2f}, temp: {:04.2f}'.format(t,1/beta))
   # draw.text((0,int(imgSize[1]/barRatio/2)),u'\u03b2: 0.4'.format(i),font=font,fill=255)
    return mplfig_to_npimage(fig)


print('starting video encoding with'+str(len(imSeq))+'frames')
animation = VideoClip(make_dframe,duration= (fN['time'].iloc[-1]))
animation.write_videofile('defect.mp4',fps = 24,audio=False,threads=4)
#danimation = VideoClip(make_dframe,duration=(len(imSeq)-1)/fps)
#danimation.write_videofile('dT.mp4',fps = 24,audio=False,threads=5)
