from __future__ import print_function
import numba
import numpy as np
import sys
import lattice_collection as lc
import tools_v1 as tool
import datetime
import params
import Relative_tools
import copy

### File with Lattice class to sweep and generate lattices and functions: 
### plaquette, average plaquette, polyakov, planar and non-planar wilson loops, wilson action,
### operator density and operator sum helper functions, topological charge included.


#-------------Measurment code -------------------
### some functions are reproduced here, outside of Lattice class, to be accessible via function call.
### - a bit redundant
def fn_periodic_link(U, txyz, direction):
    Nt, Nx, Ny, Nz = len(U), len(U[0]), len(U[0][0]), len(U[0][0][0])
    return U[txyz[0] % Nt][txyz[1] % Nx][txyz[2] %Ny][txyz[3] % Nz][direction]

def fn_move_forward_link(U, txyz, direction):
    link = fn_periodic_link(U, txyz, direction)
    new_txyz = txyz[:]
    new_txyz[direction] += 1
    return link, new_txyz

def fn_move_backward_link(U, txyz, direction):
    new_txyz = txyz[:]
    new_txyz[direction] -= 1
    link = fn_periodic_link(U, new_txyz, direction).conj().T
    return link, new_txyz

def fn_line_move_forward(U, line, txyz, direction):
    link, new_txyz = fn_move_forward_link(U, txyz, direction)
    new_line = np.dot(line, link)
    return new_line, new_txyz

def fn_line_move_backward(U, line, txyz, direction):
    link, new_txyz = fn_move_backward_link(U, txyz, direction)
    new_line = np.dot(line, link)
    return new_line, new_txyz

### plaquette calculation
def fn_plaquette(U, t, x, y, z, mu, nu):
    Nt = len(U)
    Nx = len(U[0])
    Ny = len(U[0][0])
    Nz = len(U[0][0][0])
    start_txyz = [t,x,y,z]
    result = 1.
    result, next_txyz = fn_line_move_forward(U, 1., start_txyz, mu)
    result, next_txyz = fn_line_move_forward(U, result, next_txyz, nu)
    result, next_txyz = fn_line_move_backward(U, result, next_txyz, mu)
    result, next_txyz = fn_line_move_backward(U, result, next_txyz, nu)    
    return result

### Kogut et al, PRL51 (1983) 869, Quark and gluon latent heats at the deconfinement phase transtion in SU(3) gauge theory
### energy density: \varepsilon = \beta / Nt / Ns^3 { (\sum_{space} 1 - ReTrUUUU /3 ) - (\sum{time} 1 - ReTrUUUU /3 )}
### this is just the leading term
def fn_energy_density(U, beta):
    Nt, Nx, Ny, Nz = map(len, [U, U[0], U[0][0], U[0][0][0]])    
    temporal, spatial = 0., 0.
    for t in range(Nt):
        for x in range(Nx):
            for y in range(Ny):
                for z in range(Nz):
                    for mu in range(4):
                        for nu in range(mu):
                            plaq = fn_plaquette(U, t, x, y, z, mu, nu)
                            plaq = (np.add(plaq, plaq.conj().T))       # avg both orientations averaged
                            plaq = np.trace(plaq.real) / 3. / 2.       # divide by 3 for su3 and 2 for both orientations
                            if mu == 0 or nu == 0:                     # a temporal plaquette
                                temporal += (1 - plaq) 
                            else:
                                spatial  += (1 - plaq)
    energy_dens = spatial - temporal
    energy_des = energy_dens * beta / Nt / Nx / Ny / Nz
    return energy_dens

### calculate average plaquette for u0 = <P_{\mu\nu}>^0.25
#@numba.njit
def fn_average_plaquette(U):
    Nt, Nx, Ny, Nz = map(len, [U, U[0], U[0][0], U[0][0][0]])
    res = np.zeros(np.shape(U[0,0,0,0,0,:,:]), dtype='complex128')
    for t in range(Nt):
        for x in range(Nx):
            for y in range(Ny):
                for z in range(Nz):
                    for mu in range(1, 4):
                        for nu in range(mu):
                            res = np.add(res, fn_plaquette(U, t, x, y, z, mu, nu))
    return np.trace(res).real / 3. / Nt / Nx / Ny / Nz / 6.

### Wilson action at a specific point
### S = \sum_x \sum_{\mu > \nu} (1 - 1/3 Re Tr P_{\mu\nu}(x))
###   * P_{\mu\nu}(x) = U_\mu(x) U_\nu(x + \hat\mu) U^\dagger_mu(x + \hat\nu) U^\dagger_\nu(x)
###   * fn_plaquette(U,t,x,y,z,mu,nu) returns the product of links around the plaquette, P_{\mu\nu}(x)
###   * beta = 6 / g^2
def fn_eval_point_S(U, t, x, y, z, beta, u0 = 1.):
    tmp = 0.
    for mu in range(1, 4):  #sum over \mu > \nu spacetime dimensions
        for nu in range(mu):
            tmp += ( 1. - np.real(np.trace( fn_plaquette(U, t, x, y, z, mu, nu) )) / 3. / u0**4 )
    return beta * tmp

### Calculate density for given operator.
### Requires lattice and operator to calculate along with all arguments that need to be passed to operator.
def fn_operator_density(U, operator_function, *args):
    Nt, Nx, Ny, Nz = map(len, [U, U[0], U[0][0], U[0][0][0]])
    tmp = [[[[0 for z in range(Nz)] for y in range(Ny)] for x in range(Nx)] for t in range(Nt)]
    for t in range(Nt):
        for x in range(Nx):
            for y in range(Ny):
                for z in range(Nz):
                    tmp[t][x][y][z] = operator_function(U, t, x, y, z, *args)
    return tmp


### Calculate sum for given operator over whole lattice.
### Requires lattice and operator to calculate along with all arguments that need to be passed to operator.
def fn_sum_over_lattice(U, operator_function, *args):
    Nt, Nx, Ny, Nz = map(len, [U, U[0], U[0][0], U[0][0][0]])
    sum_lattice = 0.
    for t in range(Nt):
        for x in range(Nx):
            for y in range(Ny):
                for z in range(Nz):
                    sum_lattice += operator_function(U, t, x, y, z, *args)
    return sum_lattice

### Planar Wilson loop with one dimension in time
def fn_wilson(U, t, x, y, z, mu, R, T):  #mu spatial
    #Imagine 2D loop as ABCD where A, B, C, D are edges. Split to two lines:
    #lower that consists of A*B edges
    #upper that consists of C*D edges
    #need to compute specific coordinate points where each line starts.
    if mu == 0:
        print("Wilson loop RT can not be in time direction.")
        exit()
    
    #need start of edge A and start of edge C. start of the two wilson lines
    pointA = [t, x, y, z] 
    pointC = [t, x, y, z]
    pointC[0] += T
    pointC[mu] += R
    lower = 1.
    upper = 1.
    #multiply in correct order
    for nt in range(T):
        lower, pointA = fn_line_move_forward(U, lower, pointA, 0)
        upper, pointC = fn_line_move_backward(U, upper, pointC, 0)
    for nx in range(R):
        lower, pointA = fn_line_move_forward(U, lower, pointA, mu)
        upper, pointC = fn_line_move_backward(U,upper,  pointC, mu)
    result = np.dot(lower, upper)
    return np.trace(result).real / 3.

### average of Wilson loop
def fn_wilson_average(U, R, T):
    Nt, Nx, Ny, Nz = map(len, [U, U[0], U[0][0], U[0][0][0]])
    sum_wilson = 0.
    for t in range(Nt):
        for x in range(Nx):
            for y in range(Ny):
                for z in range(Nz):
                    for direc in range(1, 4):
                        sum_wilson += fn_wilson(U, t,x,y,z, direc, R, T)
    return sum_wilson / Nx / Ny / Nz / Nt / 3.

### Wilson at a specific R.
### could be calculated with wilson operator as well
### clearer to have this option too
def fn_wilson_loops_at_r(U, R):
    Nt, Nx, Ny, Nz = map(len, [U, U[0], U[0][0], U[0][0][0]])
    wilson_loops = []
    for T in range(Nt):
        tmp = fn_wilson_average(U, R, T)
        wilson_loops.append(tmp)
    return wilson_loops

### non planar Wilson loop -> needed to get intermediate values of potential
def fn_nonplanar_wilson(U, t, x, y, z, mu, R, rho, S, T):
    #imagine non planar loop in 3D as ABCDEF where A, B, C, D, E, F are the EDGES
    #need to calculate all edges and then multiply in order A*B*C*D*E*F to get correct result
    #need to determine specific coordinate points where loop starts
    if mu == 0 or rho == 0:
        print("Non planar Wilson loop RST can not be in time direction.")
        exit()
    if mu == rho:
        print("Non planar Wilson loop RST can not have same spatial directions R, S.")
        exit()

    #Starting points for each edge. pointA for edge A etc
    pointA = [t,x,y,z]
    
    pointB = [t,x,y,z]
    pointB[0] += T
    
    pointC = [t,x,y,z]
    pointC[0] += T
    pointC[mu] += R

    pointD = [t,x,y,z]
    pointD[0] += T
    pointD[mu] += R
    pointD[rho] += S

    pointE = [t,x,y,z]
    pointE[mu] += R
    pointE[rho] += S

    pointF = [t,x,y,z]
    pointF[mu] += R

    edgeA, edgeB, edgeC, edgeD, edgeE, edgeF = 1., 1., 1., 1., 1., 1.
    for nt in range(T):
        edgeA, pointA = fn_line_move_forward(U, edgeA, pointA, 0)
        edgeD, pointD = fn_line_move_backward(U, edgeD, pointD, 0)
    
    for nr in range(R):
        edgeB, pointB = fn_line_move_forward(U, edgeB, pointB, mu)
        edgeF, pointF = fn_line_move_backward(U, edgeF, pointF,mu )

    for ns in range(S):
        edgeC, pointC = fn_line_move_forward(U, edgeC, pointC, rho)
        edgeE, pointE = fn_line_move_backward(U, edgeE, pointE,rho )

    #create full loop -> careful order required
    loop = np.dot(np.dot(np.dot(np.dot(np.dot(edgeA,
                                              edgeB),
                                              edgeC),
                                              edgeD),
                                              edgeE),
                                              edgeF)
    return np.trace(loop).real / 3.

### same as Wilson for nonplanar
def fn_nonplanar_wilson_average(U, R, S, T):
    Nt, Nx, Ny, Nz = map(len, [U, U[0], U[0][0], U[0][0][0]])
    sum_wilson = 0.
    count = 0 #for averaging
    for t in range(Nt):
        for x in range(Nx):
            for y in range(Ny):
                for z in range(Nz):
                    for mu in range(1, 4): #can not pass time in wilson R, T loop
                        for rho in range(1, 4):
                            if mu != rho:
                                sum_wilson += fn_nonplanar_wilson(U, t,x,y,z, mu, R, rho, S,  T)
                                count += 1
    return sum_wilson / float(count)

#same as wilson for nonplanar
def fn_nonplanar_wilson_loops_at_r(U, R, S):
    Nt, Nx, Ny, Nz = map(len, [U, U[0], U[0][0], U[0][0][0]])
    wilson_loops = []
    for T in range(Nt):
        tmp = fn_nonplanar_wilson_average(U, R, S, T)
        wilson_loops.append(tmp)
    return wilson_loops

### Polyakov loop 
def fn_polyakov(U):
    Nt, Nx, Ny, Nz = map(len, [U, U[0], U[0][0], U[0][0][0]])
    ans = 0.
    for x in range(Nx):
        for y in range(Ny):
            for z in range(Nz):
                p = U[0][x][y][z][0]
                for t in range(1,Nt):
                    p = np.dot(p, U[t][x][y][z][0])
                ans  = np.add(ans, np.trace(p))
    return ans / Nx / Ny / Nz 

### Polyakov density
def fn_polyakov_atpoint(U, x, y, z):
    Nt, Nx, Ny, Nz = len(U), len(U[0]), len(U[0][0]), len(U[0][0][0])
    line = 1.
    for t in range(Nt):
        line = np.dot(line, U[t][x][y][z][0])
    return np.trace(line)

### topological charge that works with only 6 terms
def fn_topological_charge(U, t, x, y, z):
    F01 = fn_F_munu(U, t, x, y, z, 0, 1)
    F23 = fn_F_munu(U, t, x, y, z, 2, 3)
    F02 = fn_F_munu(U, t, x, y, z, 0, 2)
    F31 = fn_F_munu(U, t, x, y, z, 3, 1)
    F03 = fn_F_munu(U, t, x, y, z, 0, 3)
    F12 = fn_F_munu(U, t, x, y, z, 1, 2)
    result = np.trace( np.dot(F01, F23) + np.dot(F02, F31) + np.dot(F03, F12))
    return result / ( 4. * np.pi**2 )

### antihermitian, traceless version of field strength
def fn_F_munu(U, t, x, y, z, mu, nu):
    Pmunu = fn_plaquette(U, t, x, y, z, mu, nu)
    return -1.0J * (np.subtract(Pmunu, Pmunu.conj().T) - np.trace(np.subtract(Pmunu, Pmunu.conj().T)) / 3.) / 2.


#-------------Generation code -------------------
### function called by multiprocessor in generate script
def generate(beta, u0, action, Nt, Nx, Ny, Nz, startcfg, Ncfg, thermal, border, Nhits, Nmatrix, epsilon, magnitude_1, Nu0_step='', Nu0_avg = 10):    
    
    ### loop over (t,x,y,z) and mu and set initial collection of links
    ### Either:
    ###  1. initialize to warm start by using random collection of SU(3) links, or
    ###  2. read in a previously generated configuration and continue with that Markov chain.

    name = action +'_' + str(Nt) + 'x' + str(Nx) + 'x' + str(Ny) + 'x' + str(Nz) + '_b' + str(int(beta * 100))  + '_border_' + str(border) + '_magnitude_' + str(magnitude_1)
    aa = tool.fn_a( beta ) ### returns lattice spacing calculated for the beta value
    kappa = Relative_tools.kappa_calc(aa)   ### Returns the lattice adjusted kappa value

    print('simulation parameters:')
    print('      action: ' + action)
    print('Nt,Nx,Ny,Nz = ' + str(Nt) + ',' + str(Nx) + ',' + str(Ny) + ',' + str(Nz))
    print('       beta = ' + str(beta))
    print('         u0 = ' + str(u0))
    print('      Nhits = ' + str(Nhits))
    print('      start = ' + str(startcfg))
    print('     sweeps = ' + str(Ncfg))
    print('          a = ' + str(aa) + ' fm')
    print('        1/a = ' + str(params.hbarc_GeVfm / aa) + ' GeV')
    print('        aNx = ' + str(aa * Nx) + ' fm')
    print('Temperature = ' + str(1000. * params.hbarc_GeVfm / (Nt * aa)) + ' MeV')

    if startcfg == 0:
        U = lattice(Nt, Nx, Ny, Nz, beta, u0)
    else:
        #print(action)
        U = lc.fn_load_configuration(action, Nt, Nx, Ny, Nz, beta, startcfg, border, magnitude_1, "./logs/")
        SP = lc.fn_load_spacetime(action, Nt, Nx, Ny, Nz, beta, startcfg, border, magnitude_1, "./Deformations/")
        U = lattice(Nt, Nx, Ny, Nz, beta, u0, U, SP) 
    ### I could implement something like above for spacetime deformations, but right now the code should work, and if startcfg =! 0, it recalculates the last iteration
    
    print('Continuing from cfg: ', startcfg)
    print('... generating lattices')
    matrices = Relative_tools.create_su3_set(epsilon, Nmatrix)
    acceptance = U.markov_chain_sweep(epsilon, magnitude_1, Ncfg, matrices, startcfg, name, thermal, border, Nhits, action, Nu0_step, Nu0_avg)
    print("acceptance:", acceptance)
    return acceptance



### LATTICE CLASS
class lattice():
    ### Lattice initialization.
    ### If U not passed, lattice of identities returned.
    ### Class to avoid use of incorrect initialization and passing a lot of variables
    #@numba.njit
    def __init__(self, Nt, Nx, Ny, Nz, beta, u0, U=None, SP=None):
        if None == U:
            # initialize to identities
            U = [[[[[np.identity(3, dtype='complex128') for mu in range(4)] for z in range(Nz)] for y in range(Ny)] for x in range(Nx)] for t in range(Nt)]
        # convert to numpy arrays -> significant speed up
        self.U = np.array(U)
        self.beta = beta
        self.u0 = u0
        self.Nx = Nx
        self.Ny = Ny
        self.Nz = Nz
        self.Nt = Nt
        aa = tool.fn_a( beta ) ### returns lattice spacing calculated for the beta value
        self.kappa=Relative_tools.kappa_calc(aa)
        self.aa=aa
        ### Create a separate matrix for spacetime deformations, which is the zero matrix initially (there are no deformations).
        # if None == SP:
        #     SP=np.zeros((Nt, Nx, Ny, Nz, 4), dtype='double')
        self.SP = np.array(SP, dtype='double')



        
    ### calculate link imposing periodic boundary conditions
    def periodic_link(self, txyz, direction):
        return self.U[txyz[0] % self.Nt, txyz[1] % self.Nx, txyz[2] % self.Ny, txyz[3] % self.Nz, direction, :, :]
    
    def move_forward_link(self, txyz, direction):
        link = self.periodic_link(txyz, direction)
        new_txyz = txyz[:]
        new_txyz[direction] += 1
        return link, new_txyz

    def move_backward_link(self, txyz, direction):
        new_txyz = txyz[:]
        new_txyz[direction] -= 1
        link = self.periodic_link(new_txyz, direction).conj().T
        return link, new_txyz

    def line_move_forward(self, line, txyz, direction):
        link, new_txyz = self.move_forward_link(txyz, direction)
        new_line = np.dot(line, link)
        return new_line, new_txyz

    def line_move_backward(self, line, txyz, direction):
        link, new_txyz = self.move_backward_link(txyz, direction)
        new_line = np.dot(line, link)
        return new_line, new_txyz

    ###WILSON ACTION staple
    #@numba.njit
    def dS_staple(self, t, x, y, z, mu):
        tmp = np.zeros((3, 3), dtype='complex128')
        for nu in range(4):
            if nu != mu:

                #Determine required points for the calculation of the action
                start_txyz = [t, x, y, z]
                start_txyz[mu] += 1

                ### staple 1
                line1 = 1.
                line1, next_txyz = self.line_move_forward(line1, start_txyz, nu)
                line1, next_txyz = self.line_move_backward(line1, next_txyz, mu)
                line1, next_txyz = self.line_move_backward(line1, next_txyz, nu)
                tmp += line1
                
                ### staple 2
                line2 = 1.
                line2, next_txyz = self.line_move_backward(line2, start_txyz, nu)
                line2, next_txyz = self.line_move_backward(line2, next_txyz, mu)
                line2, next_txyz = self.line_move_forward(line2, next_txyz, nu)
                tmp += line2
        
        return tmp / self.u0**3
    
    ### Improved action with rectangles
    def dS_staple_rectangle(self, t, x, y, z, mu):
        plaquette = np.zeros((3, 3), dtype = 'complex128')
        rectangle = np.zeros((3, 3), dtype = 'complex128')

        #loop through nu different than mu
        for nu in range(4):
            if nu != mu:
                start_txyz = [t, x, y, z]
                start_txyz[mu] += 1

                #positive plaquette
                line = 1.
                line, next_txyz = self.line_move_forward(line, start_txyz, nu)
                line, next_txyz = self.line_move_backward(line, next_txyz, mu)
                line, next_txyz = self.line_move_backward(line, next_txyz, nu)
                plaquette += line
                
                #negative plaquette
                line = 1.
                line, next_txyz = self.line_move_backward(line, start_txyz, nu)
                line, next_txyz = self.line_move_backward(line, next_txyz, mu)
                line, next_txyz = self.line_move_forward(line, next_txyz, nu)
                plaquette += line
                
                #rectangle Right right up left left down (Rrulld)
                #capital is the link that we compute staples around -> NOT INCLUDED IN STAPLE
                #NOTE: easier to draw individually to see what they are
                line = 1. 
                line, next_txyz = self.line_move_forward(line, start_txyz, mu)
                line, next_txyz = self.line_move_forward(line, next_txyz, nu)
                line, next_txyz = self.line_move_backward(line, next_txyz, mu)
                line, next_txyz = self.line_move_backward(line, next_txyz, mu)
                line, next_txyz = self.line_move_backward(line, next_txyz, nu)
                rectangle += line
                
                #rectangle Rulldr
                line = 1.
                line, next_txyz = self.line_move_forward(line, start_txyz, nu)
                line, next_txyz = self.line_move_backward(line, next_txyz, mu)
                line, next_txyz = self.line_move_backward(line, next_txyz, mu)
                line, next_txyz = self.line_move_backward(line, next_txyz, nu)
                line, next_txyz = self.line_move_forward(line, next_txyz, mu)
                rectangle += line

                #Ruuldd
                line = 1.
                line, next_txyz = self.line_move_forward(line, start_txyz, nu)
                line, next_txyz = self.line_move_forward(line, next_txyz, nu)
                line, next_txyz = self.line_move_backward(line, next_txyz, mu)
                line, next_txyz = self.line_move_backward(line, next_txyz, nu)
                line, next_txyz = self.line_move_backward(line, next_txyz, nu)
                rectangle += line

                #Rrdllu
                line = 1.
                line, next_txyz = self.line_move_forward(line, start_txyz, mu)
                line, next_txyz = self.line_move_backward(line, next_txyz, nu)
                line, next_txyz = self.line_move_backward(line, next_txyz, mu)
                line, next_txyz = self.line_move_backward(line, next_txyz, mu)
                line, next_txyz = self.line_move_forward(line, next_txyz, nu)
                rectangle += line
                
                #Rdllur
                line = 1.
                line, next_txyz = self.line_move_backward(line, start_txyz, nu)
                line, next_txyz = self.line_move_backward(line, next_txyz, mu)
                line, next_txyz = self.line_move_backward(line, next_txyz, mu)
                line, next_txyz = self.line_move_forward(line, next_txyz, nu)
                line, next_txyz = self.line_move_forward(line, next_txyz, mu)
                rectangle += line
                
                #Rddluu
                line = 1.
                line, next_txyz = self.line_move_backward(line, start_txyz, nu)
                line, next_txyz = self.line_move_backward(line, next_txyz, nu)
                line, next_txyz = self.line_move_backward(line, next_txyz, mu)
                line, next_txyz = self.line_move_forward(line, next_txyz, nu)
                line, next_txyz = self.line_move_forward(line, next_txyz, nu)
                rectangle += line

        ### Return staple corrected with rectangles
        return (5. * plaquette / self.u0**3 / 9.) - (rectangle / self.u0**5 / 36.)  




    ### Difference of action at a point for fixed staple. Gets link, updated link, and staple A.
    def deltaS(self, link, updated_link, staple):
        Lqcd_nu=0.
        Lqcd_old=0.
        Lqcd_nu=(self.beta * (1. - ((1. / 3.0 / self.u0) * np.real(np.trace(np.dot( (updated_link), staple))))))  ### Finds the old Lagrangian density for a point 
        Lqcd_old=(self.beta * (1. - ((1. / 3.0 / self.u0) * np.real(np.trace(np.dot( (link), staple)))))) ### Finds new Lagrangian density for a point
        ### Expanded this into two lines for no good reason, doesn't change anything, and the computational time is about the same
        ### Wanted to see whether the expanded deltaS gives back the same deltaS as the old one. It does.
        return Lqcd_nu-Lqcd_old
    
    ### Finds the Einstein - Hilbert action
    def deltaSEH(self, link, updated_link, staple, SP, SPrime, t, x, y, z, matrices, mu):
        approx_nu=Relative_tools.first_approx_tool(SPrime, t, x, y, z)      ### Returns the approximation value for the new iteration of spacetime deformations
        approx_old=Relative_tools.first_approx_tool(SP, t, x, y, z)         ### Returns the old approximation value for spacetime deformation
        ### The next three lines actually find the difference in action between old and new values of action.
        ### In hindsight, I have no idea why I've done it this way
        ### I could have saved up some lines by sending the old configuration and the new one, and finding the difference in action in the end. 
        Action_1=self.action_1(self.U[x, t, y, z, mu, :, :], updated_link, staple, approx_nu, approx_old)    ### Returns the first line of eq. 22 of spacetime note
        Action_2=self.action_2(link, updated_link, staple, self.SP, SPrime, t, x, y, z, matrices)   ### Returns the second line of eq. 22
        Action_3, Ricci=self.action_3(SP, SPrime, t, x, y, z, approx_nu, approx_old) ### Third line of eq. 22
        Action_4=self.action_4(SP, SPrime, t, x, y, z)     ### Currently returns zero for some reason. Need to investigate (might be a bug, might be a feature)
        if Action_4 != 0:
            print('Hello')
        Action=Action_1+Action_2+Action_3+Action_4       ### Finds the sum of the actions
        # print(approx_nu)
        # print(approx_old)
        # print('Action 1: ', Action_1)
        # print('Action 2: ', Action_2)
        # print('Action 3: ', Action_3)
        # print('Action 4: ', Action_4)
        # print('Ricci ', Ricci)
        # if (np.exp(-1. * Action) > np.random.uniform(0, 1)):
        #     if Ricci >= 0:
        #         print('AAAAAAAAAAAAAA', Ricci)
        # print(Action)
        return Action, Ricci

    def action_1(self, link, updated_link, staple, approx_nu, approx_old): 
        ### Expanded the difference into two lines
        ### Mostly imported from the the original code
        Lqcd_nu=0.
        Lqcd_old=0.
        Lqcd_nu=(self.beta * (1. - ((1. / 3.0 / self.u0) * np.real(np.trace(np.dot( (updated_link), staple))))))  ### Finds the old Lagrangian density for a point 
        Lqcd_old=(self.beta * (1. - ((1. / 3.0 / self.u0) * np.real(np.trace(np.dot( (link), staple)))))) ### Finds new Lagrangian density for a point
        Action_1=(1-approx_nu)*Lqcd_nu-(1-approx_old)*Lqcd_old  ### Finds the difference in the actions (first line in eq. 22)
        return Action_1
    

    ### I need to take the SPrime inside, and for each term take into account old coords (This has been resolved)
    def action_2(self, link, updated_link, staple, SP, SPrime, t, x, y, z, matrices):
        ### Returns the second line of eq. 22
        ### The fn. Plaq_approx already does everything, so this function is a bit redundant
        Plaquette_approximations=Relative_tools.Plaq_approx(self, t, x, y, z, matrices, SPrime)
        Action_2=Plaquette_approximations * self.beta
        return Action_2
    
    ### Returns third line of eq. 22.
    ### Sends the old and new spacetime configurations and finds the difference based on the equation. 
    ### Returns the updated Ricci scalar curvature as well.
    def action_3(self, SP, SPrime, t, x, y, z, approx_nu, approx_old):
        Ricci_nu=Relative_tools.Ricci(SPrime, t, x, y, z)
        Ricci_old=Relative_tools.Ricci(SP, t, x, y, z)
        Action_3=(1-approx_nu)*self.kappa*Ricci_nu-(1-approx_old)*Ricci_old*self.kappa
        # if Ricci_nu != Ricci_old:
        #     print('Hello')
        return Action_3, Ricci_nu

    ### Should return the fourth line of the Action.
    ### Given a point, finds the difference between the old Ricci value and new Ricci value in a given direction.
    ### Something's wrong
    def action_4(self, SP, SPrime, t, x, y, z):
        ###################### MULTIPLY WITH MPMATH
        # if t == 5:
        #     if x == 5:
        #         if y == 5:
        #             if z == 5:
        #                 print('Hammer time')
        coords=[t, x, y, z]
        coords_og=[t, x, y, z]
        Action_4=0.
        for alpha in range(4):
            coords[alpha]+=1
            Ricci_nu=Relative_tools.Ricci(SPrime, coords[0], coords[1], coords[2], coords[3])
            Ricci_nu_og=Relative_tools.Ricci(SPrime, coords_og[0], coords_og[1], coords_og[2], coords_og[3])
            # print('Ricci nu_og - Ricci_nu: ', Ricci_nu_og - Ricci_nu)

            Ricci_old=Relative_tools.Ricci(SP, coords[0], coords[1], coords[2], coords[3])
            Ricci_old_og=Relative_tools.Ricci(SP, coords_og[0], coords_og[1], coords_og[2], coords_og[3])
            # print('Ricci old_og - Ricci_old: ', Ricci_old_og - Ricci_old)
            # if Ricci_nu != 0:
            #     print('Ricci_nu: ', Ricci_nu)
            #     print('Ricci_old: ', Ricci_old)
            Action_4=Action_4-self.kappa*SPrime[coords[0], coords[1], coords[2], coords[3], alpha]*(Ricci_nu-Ricci_nu_og)+self.kappa*SP[coords[0], coords[1], coords[2], coords[3], alpha]*(Ricci_old-Ricci_old_og)
            coords[alpha]-=1
        return Action_4


    #@numba.njit
    def plaquette(self, t, x, y, z, mu, nu):
        Nt, Nx, Ny, Nz = self.Nt, self.Nx, self.Ny, self.Nz
        start_txyz = [t, x, y, z]
        result = 1.
        result, next_txyz = self.line_move_forward(1., start_txyz, mu)
        result, next_txyz = self.line_move_forward(result, next_txyz, nu)
        result, next_txyz = self.line_move_backward(result, next_txyz, mu)
        result, next_txyz = self.line_move_backward(result, next_txyz, nu)
        return result
    
    #@numba.njit
    def average_plaquette(self):
        Nt, Nx, Ny, Nz = self.Nt, self.Nx, self.Ny, self.Nz
        res = np.zeros(np.shape(self.U[0, 0, 0, 0, 0, :, :]), dtype='complex128')
        for t in range(Nt):
            for x in range(Nx):
                for y in range(Ny):
                    for z in range(Nz):
                        for mu in range(1, 4):
                            for nu in range(mu):
                                res = np.add(res, self.plaquette(t, x, y, z, mu, nu))
        return np.trace(res).real / 3. / Nt/ Nx / Ny / Nz / 6.
    
    
    ### Markov chain sweep. Requires: 
    ###   number of cfgs,
    ###   set of matrices to generate update,
    ###   initial cfg,
    ###   save name (if given, otherwise will not save),
    ###   hits per sweep,
    ###   action-> W for Wilson or WR for Wilson with rectangles
    ###            W_T or WR_T for tadpole improvement
    def markov_chain_sweep(self, epsilon, magnitude_1, Ncfg, matrices, initial_cfg=0, save_name='', thermal=10, border=2, Nhits=10, action='W', Nu0_step='', Nu0_avg=10):
        ratio_accept = 0.
        matrices_length = len(matrices)

        if save_name:
            output_1 = 'logs/' + save_name + '/link_' + save_name + '_'
            output_2 = 'Deformations/' + save_name + '/link_' + save_name + '_'
            output_3 = 'Rich/' + save_name + '/link_' + save_name + '_'
            ##########################
        
        #if tadpole improving, initialize list of u0 values
        if  action[-1:] == 'T':
            plaquette = []
            u0_values = [self.u0]
        
        for i in range(Ncfg - 1):
            counter_1=0
            counter_2=0
            print('starting sweep ' + str(i+initial_cfg) + ':  ' + str(datetime.datetime.now()))
            ### loop through all space time to thermalize the lattice with QCD before starting spacetime computations
            ### This part has been tested and is in agreement with the literature
            ### The spacetime part of the code simplifies to this when no spacetime perturbations are performed.
            if i < thermal and initial_cfg < thermal:   ### Checks whether the lattice has been thermalized before ~~~~##########~~~~~~~~# Somethings wrong, need to test it later on
                ### loop through spacetime dimensions
                for t in range(self.Nt):
                    for x in range(self.Nx):
                        for y in range(self.Ny):
                            for z in range(self.Nz):
                                ### loop through directions
                                for mu in range(4):
                                    ### check which staple to use
                                    if (action == 'W') or (action == 'W_T'):
                                        A =  self.dS_staple(t, x, y, z, mu) #standard Wilson or tadpole improved
                                        #(only difference is in save name of lattice since tadpole improvement is 
                                        #considered when calculating staple)
                                    elif (action == 'WR') or (action == 'WR_T'):
                                        A = self.dS_staple_rectangle(t, x, y, z, mu) #improved action with rectangles
                                        #Tadpole improve, else only half of O(a^2) error is cancelled.
                                    else:
                                        print("Error: Wrong action name or not implemented.")
                                        sys.exit()
                                    ### loop through hits
                                    for j in range( Nhits ):
                                        ### get a random SU(3) matrix
                                        r = np.random.randint(0, matrices_length)
                                        matrix = matrices[r] 
                                        ### create U'
                                        Uprime = np.dot(matrix, self.U[t, x, y, z, mu, :, :])
                                        ### calculate staple
                                        dS = self.deltaS(self.U[t, x, y, z, mu, :, :], Uprime, A)
                                        ### check if U' accepted
                                        counter_1=counter_1+1
                                        if (np.exp(-1. * dS) > np.random.uniform(0, 1)):
                                            self.U[t, x, y, z, mu, :, :] = Uprime
                                            ratio_accept += 1
                                            counter_2=counter_2 + 1
                ### Update u0. For better performance, skip every Nu0_step cfgs and append plaquettes to array. 
                ### When the array reaches size Nu0_avg, average to update u0.
                ### Wait 10 iterations from warm start.
                if action[-1:] == 'T' and (i % Nu0_step == 0) and i > 10:
                    plaquette.append( self.average_plaquette() )
                    if len(plaquette) == Nu0_avg:
                        u0_prime = ( np.mean( plaquette ) )**0.25
                        print("u0 for lattice ", self.beta, " to be updated. Previous: ", self.u0, ". New: ", u0_prime)
                        self.u0 = u0_prime           #update u0
                        u0_values.append( u0_prime ) #create array of u0 values
                        plaquette = []               #clear array

                ### save if name given
                if (save_name):
                    
                    idx = int(i) + initial_cfg
                    #print(int( idx ))
                    output_idx_1 = output_1 + str(int( idx ))
                    file_out = open(output_idx_1, 'wb')
                    np.save(file_out, self.U)  #NOTE: np.save without opening first appends .npy
                    sys.stdout.flush()
                QCD_ratio=counter_2/counter_1
                outfile = './The_ratio/' + save_name + '_QCD'
                fout = open(outfile, 'a')
                fout.write(str(i) + '   ' + str(QCD_ratio) + '\n' )
                fout.close()


        ### Same like above but with spacetime deformations and Einstein - Hilbert action
            else:
                counter_2=0
                counter_1=0
            ### loop through spacetime dimensions
                print('Lets get this ready')
                rich_array=np.zeros((self.Nt, self.Nx, self.Ny, self.Nz))   ### Initializes an array to hold the Ricci scalar curvature assigned to every point
                for t in range(self.Nt):
                    for x in range(self.Nx):
                        for y in range(self.Ny):
                            for z in range(self.Nz):
                                ### Was thinking about deforming the spacetime only once when I visit a site to potentially save up on computations
                                ### This proved to be a bit tedious to do, so decided to deform the spacetime with every hit of LQCD computations
                                ### This means that there are more spacetime deformations performed (better accuracy??)
                                ### If you can't beat them, join them
                                ### The 'for mu' part might be an issue, need to think about this

                                # Generate a spacetime grid distortion
                                ### Only add deformations if were inside the border
                                # if t >= border and t < (self.Nt-border):
                                #     if x >= border and x < (self.Nx-border):
                                #         if y >= border and y < (self.Ny-border):
                                #             if z >= border and z < (self.Nz-border):
                                                
                                ### loop through directions
                                for mu in range(4):
                                    ### check which staple to use
                                    if (action == 'W') or (action == 'W_T'):
                                        A =  self.dS_staple(t, x, y, z, mu) #standard Wilson or tadpole improved
                                        #(only difference is in save name of lattice since tadpole improvement is 
                                        #considered when calculating staple)
                                    elif (action == 'WR') or (action == 'WR_T'):
                                        A = self.dS_staple_rectangle(t, x, y, z, mu) #improved action with rectangles
                                        #Tadpole improve, else only half of O(a^2) error is cancelled.
                                    else:
                                        print("Error: Wrong action name or not implemented.")
                                        sys.exit()
                                    for j in range( Nhits ):
                                        ### get a random SU(3) matrix
                                        r = np.random.randint(0, matrices_length) 
                                        matrix = matrices[r]
                                        ### create U' and primed spacetime point
                                        Uprime = np.dot(matrix, self.U[t, x, y, z, mu, :, :])
                                        Ricci=0     ### By default, Ricci curvature of unperturbed spacetime is 0.
                                        ### The next couple of lines prepare the new spacetime (original spacetime perturbed at the given point)
                                        ### This is done only when the counter is inside the 'inner' spacetime matrix
                                        SP_prime=copy.deepcopy(self.SP[:, :, :, :, :])      ### This was needed
                                        ### Otherwise SP_prime would simply be a different name of the same memory location
                                        SP_prime=np.array(SP_prime, dtype='double')
                                        ### Need to figure out whether I need the less than or equal sign on the second parts of the conditions.
                                        ### Shouldn't matter too much anyways if the border >= 3 I think
                                        
                                        if t >= border and t < self.Nt-border:
                                            if x >= border and x < self.Nx-border:
                                                if y >= border and y < self.Ny-border:
                                                    if z >= border and z < self.Nz-border:
                                                        SP_prime[t, x, y, z, :]=SP_prime[t, x, y, z, :] + Relative_tools.Delta_gen(epsilon, magnitude_1)
                                                        dEH, Ricci = self.deltaSEH(self.U[t, x, y, z, mu, :, :], Uprime, A, self.SP, SP_prime, t, x, y, z, matrices, mu)
                                                        counter_1=counter_1+1
                                                    else:
                                                        dEH = self.deltaS(self.U[t, x, y, z, mu, :, :], Uprime, A)
                                                else:
                                                    dEH = self.deltaS(self.U[t, x, y, z, mu, :, :], Uprime, A)
                                            else:
                                                dEH = self.deltaS(self.U[t, x, y, z, mu, :, :], Uprime, A)
                                        else:
                                            dEH = self.deltaS(self.U[t, x, y, z, mu, :, :], Uprime, A)
                                        
                                        ### On the borders of the lattice, the regular action is calculated
                                        ### The same result would come out of dEH, as all the new components of the expanded action
                                        ### would get multiplied by the zero-valued space deformations

                                        ### Accepts the deformation and the updated link values only if both the change in action is accepted and the Ricci scalar curvature is positive.
                                        ### In the case when we're on the edge, still performs the test as Ricci = 0
                                        if (np.exp(-1. * dEH) > np.random.uniform(0, 1)):
                                            if Ricci >= 0:
                                               
                                                if t >= border and t < self.Nt-border:
                                                    if x >= border and x < self.Nx-border:
                                                        if y >= border and y < self.Ny-border:
                                                            if z >= border and z < self.Nz-border:
                                                                counter_2=counter_2+1
                                                self.U[t, x, y, z, mu, :, :] = Uprime
                                                self.SP[t, x, y, z, :]=copy.deepcopy(SP_prime[t, x, y, z, :])
                                                ratio_accept += 1
                                                rich_array[t, x, y, z] = Ricci
                                                
                if action[-1:] == 'T' and (i % Nu0_step == 0) and i > 10:
                    plaquette.append( self.average_plaquette() )
                    if len(plaquette) == Nu0_avg:
                        u0_prime = ( np.mean( plaquette ) )**0.25
                        print("u0 for lattice ", self.beta, " to be updated. Previous: ", self.u0, ". New: ", u0_prime)
                        self.u0 = u0_prime           #update u0
                        u0_values.append( u0_prime ) #create array of u0 values
                        plaquette = []               #clear array

                ### A bunch of saving happens below
                ### save if name given
                if (save_name):
                    idx = int(i) + initial_cfg
                    #print(int( idx ))
                    output_idx_1 = output_1 + str(int( idx ))
                    file_out = open(output_idx_1, 'wb')
                    np.save(file_out, self.U)  #NOTE: np.save without opening first appends .npy
                    sys.stdout.flush()

                    output_idx_2 = output_2 + str(int( idx ))
                    file_out_2=open(output_idx_2, 'wb')
                    np.save(file_out_2, self.SP)
                    sys.stdout.flush()

                    output_idx_3 = output_3 + str(int( idx ))
                    file_out_3=open(output_idx_3, 'wb')
                    np.save(file_out_3, rich_array)
                    sys.stdout.flush()
                
                SEH_ratio=counter_2/counter_1
                outfile = './The_ratio/' + save_name + '_SEH'
                fout = open(outfile, 'a')
                fout.write(str(i) + ' ' + str(SEH_ratio) + '\n' )

                print('Ratio of acceptence: ', SEH_ratio)
                fout.close()
        
        ############# Need to adjust the acceptance.
        ratio_accept = float(ratio_accept) / Ncfg / self.Nx / self.Ny / self.Nz / self.Nt / 4. / Nhits
        
        if action[-1:] == 'T':
            print("u0 progression: ", u0_values)
        return ratio_accept
