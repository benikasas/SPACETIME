### Script to generate quenched Wilson gauge fields, with the option of tadpole- and Symanzik-improvement.
### Modified to include multiprocessing based on Panagiotis's script.
from __future__ import print_function
import os, sys, string
import numpy as np
# from gauge_latticeqcd import *        ##I don't think its necessary
import lattice_collection as lc
from multiprocessing import Pool
import functools
import datetime
import gauge_latticeqcd

### settings
# Nt = 8
# Nx = 8
# Ny = 8
# Nz = 8
# Ncfg = 1000          # number of lattices to generate
# action = 'W'    # W = Wilson, Wilson with rectangle improvements, W_T and WR_T = With tadpole improvement
# betas = [5.5]      # betas to be generated, beta = 6/g^2
# startcfg = 900     # warm start (0) or existing cfg number to start the Markov chain
# Nhits = 50         # hits between each update
# Nmatrix = 10000    # number of random SU(3) matrices to be used for updates
# epsilon = 0.24      # how "far" away from identity the updates will be
# threads = 4       # threads used in multiprocessing

Nt = 6
Nx = 6
Ny = 6
Nz = 6
Ncfg = 1000          
action = 'W'    
betas = [5.7]      
startcfg = 0    
Nhits = 50         
Nmatrix = 10000   
epsilon = 0.20    
threads = 1 

Nu0_step = 1       # if tadpole improving, number of cfgs to skip between calculating u0.
Nu0_avg = 1        # if tadpole improving, number of u0 values to average together before updating
u0 = 1.            # u0 = <W11>^(1/4); if tadpole improving and continuing from existing lattices, set here.  Else ignore.
thermal=1         # Number of configurations before starting the general relativity part
border=2            # The border of spacetime lattice that will equal to 0

# ### initialize multiprocessing


if __name__ == '__main__':
    ### generate lattices
    for b in betas:
        dir_name = 'C:/Users/justi/SPACETIME/C_Code/logs/' + action + '_' + str(Nt) + 'x' + str(Nx) + 'x' + str(Ny) + 'x' + str(Nz) + '_b' + str(int(b * 100))
        dir_name_1 = 'C:/Users/justi/SPACETIME/C_Code/Deformations/' + action + '_' + str(Nt) + 'x' + str(Nx) + 'x' + str(Ny) + 'x' + str(Nz) + '_b' + str(int(b * 100))
        ### create output directory if it does not exist
        if not os.path.exists(dir_name):
            os.mkdir(dir_name)
        else:
            print("Directory exists for beta ", b)
        ### Same for spacetime deformations
        if not os.path.exists(dir_name_1):
            os.mkdir(dir_name_1)
        else:
            print("Deformation directory exists as well")

    #### Removed the multiprocessing as it barely does anything
    Acceptance_result=gauge_latticeqcd.generate(beta=betas[0], u0=u0, action=action, Nt=Nt, Nx=Nx, Ny=Ny, Nz=Nz, startcfg=startcfg, Ncfg=Ncfg, thermal=thermal, border=border, Nhits=Nhits, Nmatrix=Nmatrix, epsilon=epsilon, Nu0_step=Nu0_step, Nu0_avg=Nu0_avg)


    
#     p = Pool(threads)
#     # ### function to be calculated needs to use functools to work with map
#     func = functools.partial(generate, u0=u0, action=action, Nt=Nt, Nx=Nx, Ny=Ny, Nz=Nz, startcfg=startcfg, Ncfg=Ncfg, Nfluc=Nfluc, thermal=thermal, border=border, Nhits=Nhits, Nmatrix=Nmatrix, epsilon=epsilon, Nu0_step=Nu0_step, Nu0_avg=Nu0_avg)
#     p.map(func, betas) # call multiprocessing map function
#     p.terminate()      # terminate multiprocessing
