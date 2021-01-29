import sys, os, string, numba
import numpy as np
import pygal
from generate import Nt, Nx, Ny, Nz, action, epsilon
from plaq_v_cfg import Nstart, Nend
# from plaq_v_cfg import *



# data=np.loadtxt(r'plaquette_v_cfg_570_14x14x14x14_W.dat',dtype=float,delimiter=' ',skiprows=0)
# cfgno=np.loadtxt(r'plaquette_v_cfg_570_14x14x14x14_W.dat',dtype=float,delimiter=' ',skiprows=1,usecols=(0,))
# plaq=np.loadtxt(r'plaquette_v_cfg_570_14x14x14x14_W.dat',dtype=float,delimiter=' ',skiprows=1,usecols=(1,))


## Settings
# Nstart = 1
# Nend = 1790


## Uses the settings to generate the name for the selected data from relevant file.
U_infile = './dats/plaquette_from_' + str(Nstart) + '_to_' + str(Nend) + '_v_epsilon_' + str(epsilon) + '_v_cfg_' + str(int(beta * 100)) + '_' + str(Nt) + 'x' + str(Nx) + 'x' + str(Ny) + 'x' + str(Nz) + '_' + action + '.dat'

## Function to graph the plaquette values against configuration number
def grapher():
    data=np.loadtxt(U_infile,dtype=float,delimiter=' ',skiprows=0)
    rendername = './images/plaquette_from_' + str(Nstart) + '_to_' + str(Nend) + '_v_epsilon_' + str(epsilon) + '_v_cfg_' + str(int(beta * 100)) + '_' + str(Nt) + 'x' + str(Nx) + 'x' + str(Ny) + 'x' + str(Nz) + '_' + action
    lineplot=pygal.XY(stroke=False, show_legend=False, background='grey', range=(0.4, 0.6))
    lineplot.add('line', data[:])
    lineplot.render_to_file(rendername + '.svg')


## Function to run various analyses on the file
def statistician():
    plaq=np.loadtxt(U_infile,dtype=float,delimiter=' ',skiprows=1,usecols=(1,)) ## Loads plaquette data file
    totsum=np.sum(plaq) ## Finds the sum of the plaquette
    mean_val=totsum/len(plaq) ## Finds average of the plaquette
    varsum=0
    for i in range(len(plaq)):
        varsum=varsum+(plaq[i]-mean_val)**2 ## Finds the sum for the variance
    variance=1/(len(plaq)-1)*varsum ## Finds variance
    sdev=np.sqrt(variance)  ## Standard deviation of the sum
    mean_sdev=sdev/np.sqrt(len(plaq))   ## The standard deviation of the mean
    rel_sdev=sdev/mean_val  ## Relative standard deviation
    rel_error=mean_sdev/mean_val ## Relative error
    
    print('Acceptable value: ', mean_val, ' with standard deviation :', mean_sdev)


def slicer():
    in_fold='./Deformations/' + action + '_' + str(Nt) + 'x' + str(Nx) + 'x' + str(Ny) + 'x' + str(Nz) + '_b' + str(int(beta * 100))
    in_file=in_fold +  + 'link_' + action +'_' + str(Nt) + 'x' + str(Nx) + 'x' + str(Ny) + 'x' + str(Nz) + '_b' + str(int(beta * 100)) + '_' + str(100)
    SP_deformations=np.loadtxt(in_file, dtype=float)

# value=slicer()