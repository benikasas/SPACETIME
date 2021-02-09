import os, sys, string
import numba
import numpy as np
import tools_v1 as tool
import gauge_latticeqcd as gl
import params
from generate import betas, Nt, Nx, Ny, Nz, action, epsilon, u0, border, magnitude_1

### Script to calculate the evolution of the 1x1 Wilson loop as a function of Monte Carlo time



## Settings
# Nt = 8
# Nx = 8
# Ny = 8
# Nz = 8
# action = 'W'
# beta = 5.5
# u0 = 1.0
Nstart = 1
Nend = 190
beta=betas[0]

if __name__ == "__main__":

    dir = './logs/' + action + '_' + str(Nt) + 'x' + str(Nx) + 'x' + str(Ny) + 'x' + str(Nz) + '_b' + str(int(beta * 100)) + '_border_' + str(border) + '_magnitude_' + str(int(magnitude_1)) + '/'
    U_infile = dir + 'link_' + action + '_' + str(Nt) + 'x' + str(Nx) + 'x' + str(Ny) + 'x' + str(Nz) + '_b' + str(int(beta * 100)) + '_border_' + str(border) + '_magnitude_' + str(int(magnitude_1)) + '_'

    ### prepare output file
    outfile = './dats/plaquette_from_' + str(Nstart) + '_to_' + str(Nend) + '_v_epsilon_' + str(epsilon) + '_v_cfg_' + str(int(beta * 100)) + '_' + str(Nt) + 'x' + str(Nx) + 'x' + str(Ny) + 'x' + str(Nz) + '_' + action + '_border_' + str(border) + '_magnitude_' + str(int(magnitude_1)) + '.dat'
    fout = open(outfile, 'a')

    fout.write('#1:cfg  2:plaquette\n')
    for Ncfg in range(Nstart, Nend + 1):
        U = np.load(U_infile + str(Ncfg))
        pl = gl.fn_average_plaquette(U)
        fout.write(str(Ncfg) + ' ' + str(pl) + '\n' )
    #end Ncfg
    fout.close()
