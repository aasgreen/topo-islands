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
from skimage import color
import datetime
import time
import argparse
import os
import h5py
from moviepy.video.io.bindings import mplfig_to_npimage
from mpl_toolkits.mplot3d import Axes3D 
from matplotlib import cm
def read_files(args):
    '''read in, sort by time and output numpy arrays containing each image frame in the simulation'''
    
    #check if hdf5 file exists
    for fname in os.listdir(args.fName):
        if fname.endswith('h5'):
            if fname == 'hsetf.h5':
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


parser = argparse.ArgumentParser(description='Make some movies')
parser.add_argument('fName',type=str,help='directory name for data')
args = parser.parse_args()

imSeq, fN = read_files(args)
fps=24
secondsPframe_in = fN['time'].diff().mean()
fps_in = secondsPframe_in**-1
params = np.loadtxt('param.txt')
radius = np.loadtxt(args.fName+'/radius.dat')
fig = plt.figure()
ax = fig.add_subplot(111,projection = '3d')
X,Y = np.indices(imSeq[0].shape)
#a1 = ax.plot_surface(X,Y, imSeq[0], linewidth = 0)
def schler(state):
    '''take given state and return schieren state
    '''
    return np.sin(2*state)**2
def make_frame(t):
    ax.clear()
    i = int(t*fps_in)
    img = imSeq[i]
    #need to draw island
    #img = imSeq[i]
    #ax.plot_surface(X,Y,imSeq[i], linewidth = 0, cmap = cm.coolwarm, vmin = 0, vmax = 8)
    ax.imshow(img, cmap = cm.coolwarm, vmin = 0, vmax = 1)
    #ax.set_zlim([0,8])
#    a1.set_3d_properties(imSeq[i]) 
    ax.set_title(u'g: {:03.2f}, \u03b2: {:04.3f}, t: {:04.3f} '.format(params[0], params[1], t))
#    fig.canvas.draw()
    return mplfig_to_npimage(fig)

print('starting video encoding with'+str(len(imSeq))+'frames')
animation = VideoClip(make_frame,duration=(fN['time'].iloc[-1]))
animation.write_videofile('ham.mp4',fps = 24,audio=False,threads=4)
#danimation = VideoClip(make_dframe,duration=(len(imSeq)-1)/fps)
#danimation.write_videofile('dT.mp4',fps = 24,audio=False,threads=4)
