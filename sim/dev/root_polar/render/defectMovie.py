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
def read_files(args):
    '''read in, sort by time and output numpy arrays containing each image frame in the simulation'''
    
    #check if hdf5 file exists
    for fname in os.listdir(args.fName):
        if fname.endswith('h5'):
            print('h5 file found')
            hf = h5py.File(args.fName+'/'+fname)
            fN = pd.DataFrame(hf.keys())
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
    return imSeq


parser = argparse.ArgumentParser(description='Make some movies')
parser.add_argument('fName',type=str,help='directory name for data')
args = parser.parse_args()

imSeq = read_files(args)
fps=24
params = np.loadtxt('param.txt')
radius = np.loadtxt(args.fName+'/radius.dat')
def schler(state):
    '''take given state and return schieren state
    '''
    return np.sin(2*state)**2
def make_frame(t):
    i = int(t*fps)
    img = schler(imSeq[i])
    #need to draw island
    #img = imSeq[i]
    imgSize = img.shape
    xx, yy = np.indices(imgSize)
    cx, cy = [i//2 for i in imgSize]
#    print(cx,cy)
    r= np.sqrt((xx-cx)**2+(yy-cy)**2)
#    print(np.where(np.logical_and(r>=radius[i]-3, r<radius[i]+3)))
    rr=np.where(np.logical_and(r>=radius[i]-1, r<radius[i]+1))

    img[rr] = 0
#    print(radius[i])
    barRatio = 5
    #font = ImageFont.truetype("/usr/share/fonts/truetype/freefont/FreeSansBold.ttf",int(imgSize[0]/barRatio/3))
    font = ImageFont.truetype("/usr/share/fonts/truetype/freefont/FreeSansBold.ttf",int(imgSize[0]/barRatio/2))
    reSize = tuple(np.add(imgSize,(int(imgSize[0]/barRatio),0)))
#    print(reSize)
    img2 = np.ones(reSize)
    img2[int(imgSize[0]/barRatio):,:]=np.copy(img)
    im = Image.fromarray(img_as_ubyte((exposure.rescale_intensity(img2,in_range=(img.min(),img.max())))))
    draw = ImageDraw.Draw(im)
    draw.rectangle(((0,0),(imgSize[0],int(imgSize[0]/barRatio))),fill='black')
    draw.text((0,0),'g: {:03.2f}'.format(params[0]),font=font,fill=255)
    draw.text((int(imgSize[0]/2),0),u'\u03b2: {:04.3f}'.format(params[1]),font=font,fill=255)
    draw.text((0,int(imgSize[1]/barRatio/2)),u'H: {:04.3f}'.format(params[2]),font=font,fill=255)
   # draw.text((0,int(imgSize[1]/barRatio/2)),u'\u03b2: 0.4'.format(i),font=font,fill=255)
    return color.grey2rgb(np.asarray(im))

def make_dframe(t):
    i = int(t*fps)
    #img = schler(imSeq[i])
    img = dimSeq[i]
    imgSize = img.shape
    print(imgSize)
    barRatio = 5
    #font = ImageFont.truetype("/usr/share/fonts/truetype/freefont/FreeSansBold.ttf",int(imgSize[0]/barRatio/3))
    font = ImageFont.truetype("/usr/share/fonts/truetype/freefont/FreeSansBold.ttf",int(imgSize[0]/barRatio/2))
    reSize = tuple(np.add(imgSize,(int(imgSize[0]/barRatio),0)))
    print(reSize)
    img2 = np.ones(reSize)
    img2[int(imgSize[0]/barRatio):,:]=np.copy(img)
    im = Image.fromarray(img_as_ubyte((exposure.rescale_intensity(img2,in_range=(img.min(),img.max())))))
    draw = ImageDraw.Draw(im)
    draw.rectangle(((0,0),(imgSize[0],int(imgSize[0]/barRatio))),fill='black')
    draw.text((0,0),'g: {:03.2f}'.format(params[0]),font=font,fill=255)
    draw.text((int(imgSize[0]/2),0),u'\u03b2: {:04.3f}'.format(params[1]),font=font,fill=255)
    draw.text((0,int(imgSize[1]/barRatio/2)),u'H: {:04.3f}'.format(params[2]),font=font,fill=255)
   # draw.text((0,int(imgSize[1]/barRatio/2)),u'\u03b2: 0.4'.format(i),font=font,fill=255)
    return color.grey2rgb(np.asarray(im))


print('starting video encoding')
animation = VideoClip(make_frame,duration=(len(imSeq)-1)/fps)
animation.write_videofile('defect.mp4',fps = 24,audio=False,threads=4)

#danimation = VideoClip(make_dframe,duration=(len(imSeq)-1)/fps)
#danimation.write_videofile('dT.mp4',fps = 24,audio=False,threads=4)
