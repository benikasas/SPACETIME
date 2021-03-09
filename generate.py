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
# HASN'T BEEN IMPLEMENTED YET W = Wilson, Wilson with rectangle improvements, W_T and WR_T = With tadpole improvement
Nt = 11
Nx = 11
Ny = 11
Nz = 11
action = 'W'        ########## Spacetime part implemented only for Wilson action    
betas = [5.7]       # betas to be generated, beta = 6/g^2       #### For Beta=5.7, a=0.17
Nmatrix = 10000     # Number of matrices to generate for LQCD
epsilon = 0.20      # how "far" away from identity the updates will be
threads = 8         # threads used in multiprocessing
Nu0_step = 1        # if tadpole improving, number of cfgs to skip between calculating u0.
Nu0_avg = 1         # if tadpole improving, number of u0 values to average together before updating
u0 = 1.             # u0 = <W11>^(1/4); if tadpole improving and continuing from existing lattices, set here.  Else ignore.     

Ncfg = 700         # number of lattices to generate
startcfg = 299       # warm start (0) or existing cfg number to start the Markov chain
Nhits = 25          # hits between each update
thermal = 300         # Number of configurations before starting the general relativity part
border = 4            # Defines edges over which there  will be no spacetime deformations, and the action will come from LQCD only
magnitude_1 = 10**(37)  # Defines the magnitude of spacetime deformations

# ### initialize multiprocessing
if __name__ == '__main__':
    ### generate lattices
    for b in betas:
        dir_name = 'C:/Users/justi/SPACETIME/C_Code/logs/' + action + '_' + str(Nt) + 'x' + str(Nx) + 'x' + str(Ny) + 'x' + str(Nz) + '_b' + str(int(b * 100)) + '_border_' + str(border) + '_magnitude_' + str(magnitude_1)
        dir_name_1 = 'C:/Users/justi/SPACETIME/C_Code/Deformations/' + action + '_' + str(Nt) + 'x' + str(Nx) + 'x' + str(Ny) + 'x' + str(Nz) + '_b' + str(int(b * 100)) + '_border_' + str(border) + '_magnitude_' + str(magnitude_1)
        dir_name_2 = 'C:/Users/justi/SPACETIME/C_Code/Rich/' + action + '_' + str(Nt) + 'x' + str(Nx) + 'x' + str(Ny) + 'x' + str(Nz) + '_b' + str(int(b * 100)) + '_border_' + str(border) + '_magnitude_' + str(magnitude_1)
        dir_name_3 = 'C:/Users/justi/SPACETIME/C_Code/The_ratio/' + action + '_' + str(Nt) + 'x' + str(Nx) + 'x' + str(Ny) + 'x' + str(Nz) + '_b' + str(int(b * 100)) + '_border_' + str(border) + '_magnitude_' + str(magnitude_1)
        ### create output directory if it does not exist    
        if not os.path.exists(dir_name):
            os.mkdir(dir_name)
        else:
            print("Directory exists for beta ", b)
        ### Same for spacetime deformations and Ricci
        if not os.path.exists(dir_name_1):
            os.mkdir(dir_name_1)
        else:
            print("Deformation directory exists as well")
        if not os.path.exists(dir_name_2):
            os.mkdir(dir_name_2)
        else:
            print('Rich directory is rich')

    #### Removed the multiprocessing as it barely does anything
    # Acceptance_result=gauge_latticeqcd.generate(beta=betas[0], u0=u0, action=action, Nt=Nt, Nx=Nx, Ny=Ny, Nz=Nz, startcfg=startcfg, Ncfg=Ncfg, thermal=thermal, border=border, Nhits=Nhits, Nmatrix=Nmatrix, epsilon=epsilon, Nu0_step=Nu0_step, Nu0_avg=Nu0_avg)


    #### Reinstalled multiprocessing in hopes that it actually does something
    p = Pool(threads)
    ### function to be calculated needs to use functools to work with map
    func = functools.partial(gauge_latticeqcd.generate, u0=u0, action=action, Nt=Nt, Nx=Nx, Ny=Ny, Nz=Nz, startcfg=startcfg, Ncfg=Ncfg, thermal=thermal, border=border, Nhits=Nhits, Nmatrix=Nmatrix, epsilon=epsilon, magnitude_1=magnitude_1, Nu0_step=Nu0_step, Nu0_avg=Nu0_avg)
    p.map(func, betas) # call multiprocessing map function
    p.terminate()      # terminate multiprocessing
