import sys, os, string, numba
import numpy as np
import pygal
from generate import betas, Nt, Nx, Ny, Nz, action, epsilon, thermal, border, magnitude_1
from plaq_v_cfg import Nstart, Nend
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import gauge_latticeqcd
import matplotlib as cm
import tools_v1
import Relative_tools


def avrg_rich():
    end=440
    initial=-50
    dir_1 = './Rich/' + action + '_' + str(Nt) + 'x' + str(Nx) + 'x' + str(Ny) + 'x' + str(Nz) + '_b' + str(int(betas[0] * 100)) + '_border_' + str(border) + '_magnitude_' + str(int(magnitude_1)) + '/'
    U_infile = 'link_' + action + '_' + str(Nt) + 'x' + str(Nx) + 'x' + str(Ny) + 'x' + str(Nz) + '_b' + str(int(betas[0] * 100)) + '_border_' + str(border) + '_magnitude_' + str(int(magnitude_1)) + '_'
    bet=betas[0]
    aa=tools_v1.fn_a(bet)
    kappa=Relative_tools.kappa_calc(aa)
    print(kappa)
    Ricci_sum=0.
    counter=0
    for i in range(end - thermal + initial):
        idx=i+thermal-initial
        input_file_1 = dir_1 + U_infile + str(idx)
        Ricci_mat=np.load(input_file_1)
        for t in range(Nt):
            for x in range(Nx):
                for y in range(Ny):
                    for z in range(Nz):
                        if t >= border and t < Nt-border:
                            if x >= border and x < Nx-border:
                                if y >= border and y < Ny-border:
                                    if z >= border and z < Nz-border:
                                        Ricci_sum=Ricci_sum+Ricci_mat[t, x, y, z]
                                        counter=counter+1
    mean_Ricci=Ricci_sum/counter
    print('Mean Ricci scalar curvature: ', mean_Ricci)
    return mean_Ricci, counter

def cal_S(U):
    beta=betas[0]
    Nt, Nx, Ny, Nz = len(U), len(U[0]), len(U[0, 0]), len(U[0, 0, 0])
    S = 0.
    counter=0
    for t in range(Nt):
        for x in range(Nx):
            for y in range(Ny):
                for z in range(Nz):
                    if t >= border and t < Nt-border:
                        if x >= border and x < Nx-border:
                            if y >= border and y < Ny-border:
                                if z >= border and z < Nz-border:
                                    S += gauge_latticeqcd.fn_eval_point_S(U, t, x, y, z, beta, u0=1.)
                                    counter=counter+1
                #end z
            #end y
        #end x
    #end t
    return S, counter


def avrg_ac():
    Nstart = 350
    Nend = 440
    beta=betas[0]
    dir = './logs/' + action + '_' + str(Nt) + 'x' + str(Nx) + 'x' + str(Ny) + 'x' + str(Nz) + '_b' + str(int(beta * 100)) + '_border_' + str(border) + '_magnitude_' + str(int(magnitude_1)) + '/'
    U_infile = dir + 'link_' + action + '_' + str(Nt) + 'x' + str(Nx) + 'x' + str(Ny) + 'x' + str(Nz) + '_b' + str(int(beta * 100)) + '_border_' + str(border) + '_magnitude_' + str(int(magnitude_1)) + '_'
    S=0.
    counter=0
    for Ncfg in range(Nstart, Nend):
        U = np.load(U_infile + str(Ncfg))
        SS, counters=cal_S(U)
        S = S + SS
        counter=counter+counters
    mean_action=S/counter
    print('Mean QCD action: ', mean_action)
    return mean_action, counter

def sigma_ac():
    mean_action, counter=avrg_ac()
    
    Nstart = 350
    Nend = 440
    beta=betas[0]

    dir_ac = './logs/' + action + '_' + str(Nt) + 'x' + str(Nx) + 'x' + str(Ny) + 'x' + str(Nz) + '_b' + str(int(beta * 100)) + '_border_' + str(border) + '_magnitude_' + str(int(magnitude_1)) + '/'
    U_infile_ac = dir_ac + 'link_' + action + '_' + str(Nt) + 'x' + str(Nx) + 'x' + str(Ny) + 'x' + str(Nz) + '_b' + str(int(beta * 100)) + '_border_' + str(border) + '_magnitude_' + str(int(magnitude_1)) + '_'
    
    the_sum=0.
    for Ncfg in range(Nstart, Nend):
        U = np.load(U_infile_ac + str(Ncfg))
        for t in range(Nt):
            for x in range(Nx):
                for y in range(Ny):
                    for z in range(Nz):
                        if t >= border and t < Nt-border:
                            if x >= border and x < Nx-border:
                                if y >= border and y < Ny-border:
                                    if z >= border and z < Nz-border:
                                        S = gauge_latticeqcd.fn_eval_point_S(U, t, x, y, z, beta, u0=1.)
                                        the_sum=the_sum+(S-mean_action)**2
    return the_sum

def sigma_r():
    mean_Ricci, counter=avrg_rich()
    Nstart = 350
    Nend = 440
    beta=betas[0]
    dir_ric = './Rich/' + action + '_' + str(Nt) + 'x' + str(Nx) + 'x' + str(Ny) + 'x' + str(Nz) + '_b' + str(int(betas[0] * 100)) + '_border_' + str(border) + '_magnitude_' + str(int(magnitude_1)) + '/'
    U_infile_ric = dir_ric + 'link_' + action + '_' + str(Nt) + 'x' + str(Nx) + 'x' + str(Ny) + 'x' + str(Nz) + '_b' + str(int(betas[0] * 100)) + '_border_' + str(border) + '_magnitude_' + str(int(magnitude_1)) + '_'
    the_sum=0.
    for Ncfg in range(Nstart, Nend):
        Ricci_mat=np.load(U_infile_ric + str(Ncfg))
        for t in range(Nt):
            for x in range(Nx):
                for y in range(Ny):
                    for z in range(Nz):
                        if t >= border and t < Nt-border:
                            if x >= border and x < Nx-border:
                                if y >= border and y < Ny-border:
                                    if z >= border and z < Nz-border:
                                        Ricci = Ricci_mat[t, x, y, z]
                                        the_sum=the_sum+(Ricci-mean_Ricci)**2
    return the_sum




def corr_coff():
    mean_action, counter_1=avrg_ac()
    mean_Ricci, counter_2=avrg_rich()
    sigmasigma=np.sqrt(sigma_ac()*sigma_r())
    Nstart = 350
    Nend = 440
    beta=betas[0]
    counter=0

    dir_ac = './logs/' + action + '_' + str(Nt) + 'x' + str(Nx) + 'x' + str(Ny) + 'x' + str(Nz) + '_b' + str(int(beta * 100)) + '_border_' + str(border) + '_magnitude_' + str(int(magnitude_1)) + '/'
    U_infile_ac = dir_ac + 'link_' + action + '_' + str(Nt) + 'x' + str(Nx) + 'x' + str(Ny) + 'x' + str(Nz) + '_b' + str(int(beta * 100)) + '_border_' + str(border) + '_magnitude_' + str(int(magnitude_1)) + '_'
    dir_ric = './Rich/' + action + '_' + str(Nt) + 'x' + str(Nx) + 'x' + str(Ny) + 'x' + str(Nz) + '_b' + str(int(betas[0] * 100)) + '_border_' + str(border) + '_magnitude_' + str(int(magnitude_1)) + '/'
    U_infile_ric = dir_ric + 'link_' + action + '_' + str(Nt) + 'x' + str(Nx) + 'x' + str(Ny) + 'x' + str(Nz) + '_b' + str(int(betas[0] * 100)) + '_border_' + str(border) + '_magnitude_' + str(int(magnitude_1)) + '_'
    the_sum=0.
    for Ncfg in range(Nstart, Nend):
        U = np.load(U_infile_ac + str(Ncfg))
        Ricci_mat=np.load(U_infile_ric + str(Ncfg))
        for t in range(Nt):
            for x in range(Nx):
                for y in range(Ny):
                    for z in range(Nz):
                        if t >= border and t < Nt-border:
                            if x >= border and x < Nx-border:
                                if y >= border and y < Ny-border:
                                    if z >= border and z < Nz-border:
                                        S = gauge_latticeqcd.fn_eval_point_S(U, t, x, y, z, beta, u0=1.)
                                        Ricci = Ricci_mat[t, x, y, z]
                                        the_sum=the_sum+(S-mean_action)*(Ricci-mean_Ricci)
                                        counter=counter+1
    the_coff=the_sum/sigmasigma
    stdev=(1-the_coff**2)/counter
    print('The correlation coefficient: ', the_coff, ' with the standard error: ', stdev)
    print(counter)
    return the_coff



value=corr_coff()