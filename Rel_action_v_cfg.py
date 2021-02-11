import os, sys, string
import numba
import numpy as np
import tools_v1 as tool
import gauge_latticeqcd as gl
import params
from generate import Nx, Ny, Nz, Nt, action, u0, betas, border, magnitude_1, thermal


### Script to calculate the evolution of the action as a function of Monte Carlo time
Nstart = 0
Nend = 60

beta=betas[0]
#@numba.njit
def calc_S_og(U):
    Nt, Nx, Ny, Nz = len(U), len(U[0]), len(U[0, 0]), len(U[0, 0, 0])
    S = 0.
    for t in range(Nt):
        for x in range(Nx):
            for y in range(Ny):
                for z in range(Nz):

                    S += gl.fn_eval_point_S(U, t, x, y, z, beta, u0=1.)
    return S

def calc_S_rel(U, SP):
    Nt, Nx, Ny, Nz = len(U), len(U[0]), len(U[0, 0]), len(U[0, 0, 0])
    S = 0.
    for t in range(Nt):
        for x in range(Nx):
            for y in range(Ny):
                for z in range(Nz):
                    coords=[t, x, y, z]
                    S_og = gl.fn_eval_point_S(U, t, x, y, z, beta, u0=1.)
                    S_rel = 0
                    for alpha in range(4):
                        coords[alpha]+=1
                        S_rel = S_rel + (rel_fn_eval_point_S(U, coords, beta, u0=1) - S_og) * SP[t, x, y, z, alpha]
                        coords[alpha]-=1
                    S = S + S_og - S_rel
    return S

def rel_fn_eval_point_S(U, coords, beta, u0 = 1.):
    tmp = 0.
    t = coords[0]
    x = coords[1]
    y = coords[2]
    z = coords[3]
    for mu in range(1, 4):  #sum over \mu > \nu spacetime dimensions
        for nu in range(mu):
            tmp += ( 1. - np.real(np.trace( gl.fn_plaquette(U, t, x, y, z, mu, nu) )) / 3. / u0**4 )
    return beta * tmp

dir_1 = './logs/' + action + '_' + str(Nt) + 'x' + str(Nx) + 'x' + str(Ny) + 'x' + str(Nz) + '_b' + str(int(beta * 100)) + '_border_' + str(border) + '_magnitude_' + str(int(magnitude_1)) + '/'
dir_2 = './Deformations/' + action + '_' + str(Nt) + 'x' + str(Nx) + 'x' + str(Ny) + 'x' + str(Nz) + '_b' + str(int(beta * 100)) + '_border_' + str(border) + '_magnitude_' + str(int(magnitude_1)) + '/'

U_infile_1 = dir_1 + 'link_' + action + '_' + str(Nt) + 'x' + str(Nx) + 'x' + str(Ny) + 'x' + str(Nz) + '_b' + str(int(beta * 100)) + '_border_' + str(border) + '_magnitude_' + str(int(magnitude_1)) + '_'
U_infile_2 = dir_2 + 'link_' + action + '_' + str(Nt) + 'x' + str(Nx) + 'x' + str(Ny) + 'x' + str(Nz) + '_b' + str(int(beta * 100)) + '_border_' + str(border) + '_magnitude_' + str(int(magnitude_1)) + '_'


### prepare output file
outfile = './dats/rel_S_v_cfg_' + str(int(beta * 100)) + '_' + str(Nt) + 'x' + str(Nx) + 'x' + str(Ny) + 'x' + str(Nz) + '_' + action + '_border_' + str(border) + '_magnitude_' + str(int(magnitude_1)) + '.dat'
fout = open(outfile, 'a')

fout.write('#1:cfg      2:S     3:Rel_S\n')
for Ncfg in range(Nstart, Nend + 1):
    S_2=0
    U = np.load(U_infile_1 + str(Ncfg))
    S_1 = calc_S_og(U)
    if Ncfg >= thermal:
        SP = np.load(U_infile_2 + str(Ncfg))
        S_2 = calc_S_rel(U, SP)
    fout.write(str(Ncfg) + ' ' + str(S_1) + ' ' + str(S_2) + '\n' )
#end Ncfg
fout.close()
