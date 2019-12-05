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
from scipy.interpolate import griddata

def binMatrix(grid,step):
    x,y = grid.shape
    returnGrid = np.zeros([x//step,y//step])
    for i in np.arange(0,x//step):
        for j in np.arange(0,y//step):
    #        print(i,j)
            returnGrid[i,j]=grid[ step*i:step*i+step, step*j:step*j+step].mean()
    return returnGrid

def xPos(theta):
    return np.cos(theta)
def yPos(theta):
    return np.sin(theta)

def hamxy(grid):
    imShape = grid.shape
    energy = np.zeros(imShape)
    for i,row in enumerate(grid):
        for j,col in enumerate(row):
            #print(i,j)
            energy[i,j]= np.cos(col-grid[(i+1)%imShape[0],j]) + np.cos(col-grid[(i-1)%imShape[0],j])+np.cos(col-grid[i,(j+1)%imShape[0]])+np.cos(col-grid[i,(j+1)%imShape[0]])
    return -energy

def wNum(theta1,theta2):
    delA = theta1-theta2
    return delA-np.arcsin(np.sin(theta1-theta2))

def dCount(grid):
    imShape = grid.shape
    N = imShape[0]
    dgrid = np.zeros(imShape)
    wgrid = np.zeros(imShape)
    for i,row in enumerate(grid):
        for j,col in enumerate(row):
            wN = wNum(grid[i, (j-1)%N],col) +\
                 wNum(grid[ (i-1)%N, (j-1)%N], grid[i,(j-1)%N])+\
                 wNum(grid[ (i-1)%N, (j)%N],grid[ (i-1)%N, (j-1)%N])+\
                 wNum(col,grid[(i-1)%N, (j)%N])
            wgrid[i,j] = wN
            if (wN >= np.pi):
                dgrid[i,j] = 1
            elif (-wN >= np.pi):
                dgrid[i,j] = -1
    return (wgrid, dgrid)

def schler(grid):
    return np.sin(2*np.mod(grid,2*np.pi))**2
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Make some movies')
    parser.add_argument('fName',type=str,help='directory name for data')
    args = parser.parse_args()
    fileNames = glob.glob(args.fName+'/out*.dat')
    dfileNames = glob.glob(args.fName+'/defect*.dat')
#j    params = np.loadtxt('param.txt')
    fN = pd.DataFrame(fileNames)
    fN['time'] = fN[0].str.extract(r'(\d*\.?\d*)\.dat').astype('float')
    fN.sort_values(by=['time'],inplace=True)

    dfN = pd.DataFrame(dfileNames)
    dfN['time'] = dfN[0].str.extract(r'(\d*\.?\d*)\.dat').astype('float')
    dfN.sort_values(by=['time'],inplace=True)
    tPattern = re.compile('.*beta-(\d*\.?\d*)')
    beta = float(tPattern.search(args.fName)[1])
    #Sort fileNames by number

    imSeq = [np.loadtxt(f) for f in fN[0]]
    dimSeq = [np.loadtxt(f) for f in dfN[0]]
    
    f1 = binMatrix(imSeq[0],2)
    imShape = f1.shape
    X,Y = np.indices(imShape)

    fig,[ax,ax2] = plt.subplots(nrows=1,ncols=2)
    ax.set_aspect(aspect='equal')
    ax.axis('off')
    ax2.set_aspect(aspect='equal')
    ax2.axis('off')

    #ax.quiver(X,Y, xPos(imSeq[0]), yPos(imSeq[0]))

    circ = {}
    for i,j in zip(X.ravel(),Y.ravel()):
        print(i,j)
        circ[(i,j)]=(plt.Circle((i,j),.1,zorder=2))
        ax.add_artist(circ[(i,j)])

    dCirc = {}
    dText = {}
    for i, j in zip(X[:, :].ravel(), Y[:, :].ravel()):
        print(i, j)
        dCirc[(i,j)]=(plt.Rectangle(
            (i-1, j-1), 1, 1, alpha=.1, edgecolor='k', facecolor='white', zorder=0) )
        ax.add_artist(dCirc[(i,j)])
        dText[(i,j)]=(plt.Text(x=i-.5, y=j-.5, text='',fontsize=12, va='center', ha='center', zorder=1))
        ax.add_artist(dText[(i,j)])

    #f1 = binMatrix(imSeq[0],10)
    fMain = imSeq[0]
    U = xPos(f1)
    V = yPos(f1)
    gg = hamxy(f1)
    wl,dL = dCount(f1)
    programD = dimSeq[0] 

    pD = np.where(dL.ravel()==1)[0]
    mD = np.where(dL.ravel()==-1)[0]


    xpDefect = X.ravel()[pD]
    ypDefect = Y.ravel()[pD]
    xmDefect = X.ravel()[mD]
    ymDefect = Y.ravel()[mD]

    pDefect = np.vstack([xpDefect,ypDefect]).T
    mDefect = np.vstack([xmDefect,ymDefect]).T

    #find defects from the program

    progpD = np.where(programD.ravel()==1)[0]
    progmD = np.where(programD.ravel()==-1)[0]


    pxpDefect = X.ravel()[progpD]
    pypDefect = Y.ravel()[progpD]
    pxmDefect = X.ravel()[progmD]
    pymDefect = Y.ravel()[progmD]

    progpDefect = np.vstack([pxpDefect,pypDefect]).T
    progmDefect = np.vstack([pxmDefect,pymDefect]).T

    [dCirc[tuple(loc)].set_facecolor('blue') for loc in mDefect]
    [dCirc[tuple(loc)].set_facecolor('red') for loc in pDefect]


    [dText[tuple(loc)].set_text('+') for loc in progpDefect]
    [dText[tuple(loc)].set_text('-') for loc in progmDefect]

        
    Q=ax.quiver(X,Y, U,V,gg,scale=30,cmap='copper')
    #maxPos = X.max()

    #sgrid_x,sgrid_y = np.mgrid[0:maxPos:200j,0:maxPos:200j]
    #points = np.vstack([X.ravel(),Y.ravel()]).T
    #Jsgrid = griddata(points, f1.ravel(), (sgrid_x,sgrid_y),method='nearest')
    ax2.imshow(schler(fMain).T,cmap='gray',origin='lower')


    #choose frame jump
    def get_frame(i):
        return (i)%len(imSeq)
    norm = mpl.colors.Normalize(vmin=gg.min(),vmax=gg.max())
    cmap = mpl.cm.copper
    m = mpl.cm.ScalarMappable(norm=norm,cmap=cmap)
    def update_quiver(num,Q,X,Y):

        frame = imSeq[get_frame(num)]
        fMain = frame
        frame = binMatrix(frame,2)
        dframe = dimSeq[get_frame(num)]
        U = xPos(frame)
        V = yPos(frame)
        en = hamxy(frame)
        Q.set_UVC(U,V,en)

        #now, set defects

        #first, set everything back to white
        [dCirc[tuple(loc)].set_facecolor('white') for loc in list(dCirc.keys())]


        [dText[tuple(loc)].set_text('') for loc in list(dText.keys())]

        
        w, d = dCount(frame)

        pD = np.where(d.ravel() == 1)[0]
        mD = np.where(d.ravel() == -1)[0]

        
        xpDefect = X.ravel()[pD]
        ypDefect = Y.ravel()[pD]
        xmDefect = X.ravel()[mD]
        ymDefect = Y.ravel()[mD]

        pDefect = np.vstack([xpDefect,ypDefect]).T
        mDefect = np.vstack([xmDefect,ymDefect]).T


        #read in defects from program tagging
        programD =dframe 

        progpD = np.where(programD.ravel()==1)[0]
        progmD = np.where(programD.ravel()==-1)[0]


        pxpDefect = X.ravel()[progpD]
        pypDefect = Y.ravel()[progpD]
        pxmDefect = X.ravel()[progmD]
        pymDefect = Y.ravel()[progmD]

        progpDefect = np.vstack([pxpDefect,pypDefect]).T
        progmDefect = np.vstack([pxmDefect,pymDefect]).T


        [dCirc[tuple(loc)].set_facecolor('blue') for loc in mDefect]
        [dCirc[tuple(loc)].set_facecolor('red') for loc in pDefect]


        [dText[tuple(loc)].set_text('+') for loc in progpDefect]
        [dText[tuple(loc)].set_text('-') for loc in progmDefect]

        
        ax.set_title('t: {}, temp: {}'.format(get_frame(num),1/beta))

        ax2.clear()

       # sgrid_x,sgrid_y = np.mgrid[0:maxPos:200j,0:maxPos:200j]
       #J points = np.vstack([X.ravel(),Y.ravel()]).T
        #sgrid = np.mod(griddata(points, frame.ravel(), (sgrid_x, sgrid_y), method='nearest' ), 2*np.pi)

        ax2.imshow(schler(fMain).T,cmap='gray',origin='lower')

        return Q,

    anim = animation.FuncAnimation(fig,update_quiver,fargs=(Q,X,Y),interval=300,blit=False)

    fig.tight_layout

    anim.save('b-{}monte-carlo-xy-texture.avi'.format(beta))
