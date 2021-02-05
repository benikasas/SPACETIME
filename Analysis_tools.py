import sys, os, string, numba
import numpy as np
import pygal
from generate import betas, Nt, Nx, Ny, Nz, action, epsilon, thermal, border
from plaq_v_cfg import Nstart, Nend
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import gauge_latticeqcd
import matplotlib as cm
# from plaq_v_cfg import *


beta=betas[0]
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

### Used to reduce 4D to 3D, and use color grading to find the magnitude of deformations

def SP_slicer():
    end=230
    beta=5.7
    dir_1 = './Rich/' + action + '_' + str(Nt) + 'x' + str(Nx) + 'x' + str(Ny) + 'x' + str(Nz) + '_b' + str(int(betas[0] * 100)) + '/'
    dir_2 = './Rich_images/' + action + '_' + str(Nt) + 'x' + str(Nx) + 'x' + str(Ny) + 'x' + str(Nz) + '_b' + str(int(betas[0] * 100)) + '/'
    U_infile = 'link_' + action + '_' + str(Nt) + 'x' + str(Nx) + 'x' + str(Ny) + 'x' + str(Nz) + '_b' + str(int(betas[0] * 100)) + '_'
    for i in range(end - thermal):

        idx=i+thermal
        input_file = dir_1 + U_infile + str(idx)
        SP=np.load(input_file)
        x=range(len(SP))
        y=range(len(SP))
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        X, Y = np.meshgrid(x, y) 
        ax.plot_wireframe(X, Y, SP[border+1][:][:][border+1])
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        # plt.show()
        if not os.path.exists(dir_2):
            os.mkdir(dir_2)
        out_file= dir_2 + U_infile + '_' + str(idx)
        fig.savefig(out_file)

def U_slicer():
    end=230
    dir_1 = './logs/' + action + '_' + str(Nt) + 'x' + str(Nx) + 'x' + str(Ny) + 'x' + str(Nz) + '_b' + str(int(betas[0] * 100)) + '/'
    dir_2 = './Plaq_images/' + action + '_' + str(Nt) + 'x' + str(Nx) + 'x' + str(Ny) + 'x' + str(Nz) + '_b' + str(int(betas[0] * 100)) + '/'
    U_infile = 'link_' + action + '_' + str(Nt) + 'x' + str(Nx) + 'x' + str(Ny) + 'x' + str(Nz) + '_b' + str(int(betas[0] * 100)) + '_'
    for i in range(end - thermal):

        idx=i+thermal
        input_file = dir_1 + U_infile + str(idx)
        U=np.load(input_file)
        x=range(len(U))
        y=range(len(U))
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        X, Y = np.meshgrid(x, y)
        grid=np.zeros((len(U),len(U))) 
        for alpha in range(len(U)):
            for beta in range(len(U)):
                grid[alpha][beta]=gauge_latticeqcd.fn_eval_point_S(U, border+1, alpha, beta, border+1, betas[0], u0=1.)

        # norm = plt.Normalize(grid.min(), grid.max())
        # print(norm)
        # colors = plt.viridis()
        # rcount, ccount, _ = colors.shape
        
        ax.plot_surface(X, Y, grid, cmap=cm.coolwarm, linewidth=0, antialiased=False)
        ax.set_facecolor((0,0,0,0))
        plt.show()
        # if not os.path.exists(dir_2):
        #     os.mkdir(dir_2)
        # out_file= dir_2 + U_infile + '_' + str(idx)
        # fig.savefig(out_file)

def U_SP_fusion():
    end=230
    dir_1 = './logs/' + action + '_' + str(Nt) + 'x' + str(Nx) + 'x' + str(Ny) + 'x' + str(Nz) + '_b' + str(int(betas[0] * 100)) + '/'
    dir_2 = './Rich/' + action + '_' + str(Nt) + 'x' + str(Nx) + 'x' + str(Ny) + 'x' + str(Nz) + '_b' + str(int(betas[0] * 100)) + '/'
    dir_3 = './Fusion/' + action + '_' + str(Nt) + 'x' + str(Nx) + 'x' + str(Ny) + 'x' + str(Nz) + '_b' + str(int(betas[0] * 100)) + '/'
    U_infile = 'link_' + action + '_' + str(Nt) + 'x' + str(Nx) + 'x' + str(Ny) + 'x' + str(Nz) + '_b' + str(int(betas[0] * 100)) + '_'
    for i in range(end - thermal):
        idx=i+thermal
        input_file_1 = dir_1 + U_infile + str(idx)
        U=np.load(input_file_1)

        input_file_2 = dir_2 + U_infile + str(idx)
        SP=np.load(input_file_2)

        x=range(len(U))
        y=range(len(U))
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        X, Y = np.meshgrid(x, y)
        grid=np.zeros((len(U),len(U))) 
        for alpha in range(len(U)):
            for beta in range(len(U)):
                grid[alpha][beta]=gauge_latticeqcd.fn_eval_point_S(U, border+1, alpha, beta, border+1, betas[0], u0=1.)

        ax.plot_wireframe(X, Y, grid)
        ax.plot_wireframe(X, Y, 16+SP[border+1][:][:][border+1]*10**16, color='red')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        plt.show()
        if not os.path.exists(dir_3):
            os.mkdir(dir_3)
        out_file= dir_3 + U_infile + '_' + str(idx)
        fig.savefig(out_file)



value=U_slicer()

# in_fold='./Deformations/' + action + '_' + str(Nt) + 'x' + str(Nx) + 'x' + str(Ny) + 'x' + str(Nz) + '_b' + str(int(beta * 100))
# in_file=in_fold +  + 'link_' + action +'_' + str(Nt) + 'x' + str(Nx) + 'x' + str(Ny) + 'x' + str(Nz) + '_b' + str(int(beta * 100)) + '_' + str(100)
# SP_deformations=np.loadtxt(in_file, dtype=float)