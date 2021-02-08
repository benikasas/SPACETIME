import numpy as np
#from generate import Nt, Nx, Ny, Nz, action, epsilon, u0
import tools_v1



### The SU2 and SU3 computations were transported here from gauge_latticeqcd.py, from the given code

### Generate SU(2) matrix as described in Gattringer & Lang
def matrix_su2(epsilon = 0.2):
    ### Pauli matrices
    sigma1 = np.array([[0, 1], [1, 0]])
    sigma2 = np.array([[0, -1J], [1J, 0]])
    sigma3 = np.array([[1, 0], [0, -1]])
    r = [0., 0., 0., 0.]
    for i in range(4):
        r[i] = (np.random.uniform(0, 0.5))
    ### normalize
    norm = np.sqrt(r[1]**2 + r[2]**2 + r[3]**2)
    r[1:] = map(lambda x: epsilon*x / norm, r[1:])
    r[0]  = np.sign(r[0]) * np.sqrt(1. - epsilon**2)
    M = np.identity(2, dtype='complex128')
    M = M * r[0]
    M = np.add(1J * r[1] * sigma1, M)
    M = np.add(1J * r[2] * sigma2, M)
    M = np.add(1J * r[3] * sigma3, M)
    return M

### Use SU(2) matrices to generate SU(3) matrix
### From Gattringer & Lang's textbook.
### Need 3 SU(2) matrices for one SU(3) matrix
def matrix_su3(epsilon = 0.2):
    R_su2 = matrix_su2(epsilon)
    S_su2 = matrix_su2(epsilon)
    T_su2 = matrix_su2(epsilon)
    # initialise to identity, need complex numbers from now
    R = np.identity(3, dtype='complex128')
    S = np.identity(3, dtype='complex128')
    T = np.identity(3, dtype='complex128')
    # upper
    R[:2,:2] = R_su2
    # edges
    S[0:3:2, 0:3:2] = S_su2
    # lower
    T[1:,1:] = T_su2
    # create final matrix
    X = np.dot(R, S)
    return np.dot(X, T) 

### Create set of SU(3) matrices
### Needs to be large enough to cover SU(3)
def create_su3_set(epsilon = 0.2, tot = 1000):
    matrices = []
    for i in range(tot):
        X = matrix_su3(epsilon)
        matrices.append(X)
        matrices.append(X.conj().T)
    return matrices

### Dimensionless Kappa calculator
def kappa_calc(aa):
    kappa=7.6*10**37*aa*aa  ### This was found using the kappa value (from constants) multiplied by lattice spacing^2
    return kappa

## A function to generate a single disturbance in Spacetime
## Need to adjust the magnitude of epsilon
def Delta_gen(epsilon, magnitude_1):
    Delta=[0., 0., 0., 0.,]
    magnitude_2=epsilon/(magnitude_1)     ### Specify the magnitude of deformations
    for i in range(len(Delta)): 
        Delta[i]=np.random.uniform(-magnitude_2, magnitude_2) 
    return Delta

## First order approximations of spacetime deformations in all directions
## Make a for loop in the future
## Could be combined with inv_jack
## This finds the difference between the deformation at x' and x'- a * nu where nu is the direction vector
def first_approx_tool(SP, t, x, y, z):
    SP_og=SP[t, x, y, z, :]
    SP_t=SP[t+1, x, y, z, :]
    SP_x=SP[t, x+1, y, z, :]
    SP_y=SP[t, x, y+1, z, :]
    SP_z=SP[t, x, y, z+1, :]        
    diff_t=SP_t-SP_og
    diff_x=SP_x-SP_og
    diff_y=SP_y-SP_og
    diff_z=SP_z-SP_og     
    total_diff=diff_t[0]+diff_x[1]+diff_y[2]+diff_z[3]  ### Esentially finds the trace (this is the sum of delta_ro Delta^ro)
    return total_diff

### First order approximation of the Jacobian (using equation 17)
### Esentially its a deformation matrix with an identity added
### Could add a for loop as well
def inv_Jack(SPrime, t, x, y, z):
    SP_og=SPrime[t, x, y, z, :]
    SP_t=SPrime[t+1, x, y, z, :]
    SP_x=SPrime[t, x+1, y, z, :]
    SP_y=SPrime[t, x, y+1, z, :]
    SP_z=SPrime[t, x, y, z+1, :]
    jack=np.identity(4)
    jack[0][:]=np.add(jack[0][:], (SP_t-SP_og))
    jack[1][:]=np.add(jack[1][:], (SP_x-SP_og))
    jack[2][:]=np.add(jack[2][:], (SP_y-SP_og))
    jack[3][:]=np.add(jack[3][:], (SP_z-SP_og))
    jack=np.transpose(jack)     ### I am pretty sure this is needed, but I might be wrong
    jack=np.linalg.inv(jack)
    return jack

### Finds the h matrix for the S_Prime matrix with given coordinates
### This is found using equation 23
def h_matrix_producinator(SPrime, coords):
    t=coords[0]
    x=coords[1]
    y=coords[2]
    z=coords[3]
    jack=inv_Jack(SPrime, t, x, y, z)
    h_matrix=np.zeros((4, 4))
    ### Diagonal terms
    h_matrix[0, 0]=jack[0][0]*jack[0][0]+jack[1][0]*jack[1][0]+jack[2][0]*jack[2][0]+jack[3][0]*jack[3][0]-1
    h_matrix[1, 1]=jack[0][1]*jack[0][1]+jack[1][1]*jack[1][1]+jack[2][1]*jack[2][1]+jack[3][1]*jack[3][1]-1
    h_matrix[2, 2]=jack[0][2]*jack[0][2]+jack[1][2]*jack[1][2]+jack[2][2]*jack[2][2]+jack[3][2]*jack[3][2]-1
    h_matrix[3, 3]=jack[0][3]*jack[0][3]+jack[1][3]*jack[1][3]+jack[2][3]*jack[2][3]+jack[3][3]*jack[3][3]-1

    ### I am pretty sure the matrix is going to be symmetric.
    h_matrix[0, 1]=jack[0][0]*jack[0][1]+jack[1][0]*jack[1][1]+jack[2][1]*jack[2][1]+jack[3][0]*jack[3][1]
    h_matrix[0, 2]=jack[0][0]*jack[0][2]+jack[1][0]*jack[1][2]+jack[2][2]*jack[2][2]+jack[3][0]*jack[3][2]
    h_matrix[0, 3]=jack[0][0]*jack[0][3]+jack[1][0]*jack[1][3]+jack[2][3]*jack[2][3]+jack[3][0]*jack[3][3]
    h_matrix[1, 0]=h_matrix[0, 1]
    h_matrix[2, 0]=h_matrix[0, 2]
    h_matrix[3, 0]=h_matrix[0, 3]


    h_matrix[1, 2]=jack[0][1]*jack[0][2]+jack[1][1]*jack[1][2]+jack[2][1]*jack[2][2]+jack[3][1]*jack[3][2]
    h_matrix[1, 3]=jack[0][1]*jack[0][3]+jack[1][1]*jack[1][3]+jack[2][1]*jack[2][3]+jack[3][1]*jack[3][3]
    h_matrix[2, 1]=h_matrix[1, 2]
    h_matrix[3, 1]=h_matrix[1, 3]

    h_matrix[2, 3]=jack[0][2]*jack[0][3]+jack[1][2]*jack[1][3]+jack[2][2]*jack[2][3]+jack[3][2]*jack[3][3]
    h_matrix[3, 2]=h_matrix[2, 3]
    return h_matrix

### Collects surrounding h matrices for given coordinates into a 4x4 matrix 
def h_matrix_collectinator(SPrime, t, x, y, z):
    h_alphabeta=[[np.identity(4) for alph in range(4)] for beth in range(4)]    ### For h matrices where we add two lattice spacings from the origin
    h_alpha=[np.identity(4) for alph in range(4)]   ### For h matrices separated from origin by single lattice spacing
    coords=[t, x, y, z]
    for alpha in range(4):
        for beta in range(4):
            coords[alpha]=+1
            coords[beta]+=1
            h_alphabeta[alpha][beta]=h_matrix_producinator(SPrime, coords)
            coords[beta]-=1
        h_alpha[alpha]=h_matrix_producinator(SPrime, coords)
        coords[alpha]-=1
    h_og=h_matrix_producinator(SPrime, coords)      ### The h matrix for origin
    return h_alphabeta, h_alpha, h_og


### Function to find the Ricci scalar curvature at a given point
def Ricci(SPrime, t, x, y, z):
    h_alphabeta, h_alpha, h_og = h_matrix_collectinator(SPrime, t, x, y, z)
    h_alphabeta=np.array(h_alphabeta)
    h_alpha=np.array(h_alpha)
    h_og=np.array(h_og)
    R_alphabeta=0
    R_alpha=0
    R_og=0
    R_2_betabeta=0
    R_betabeta=0
    R_og_betabeta=0
    Ricci_scalar=0
    # ### First term of eq 26
    # for alpha in range(4):
    #     for beta in range(4):
    #         for i in range(4):
    #             for j in range(4):
    #                 R_alphabeta=R_alphabeta+h_alphabeta[alpha, beta, i, j]
    for alpha in range(4):
        for beta in range(4):
            R_alphabeta=R_alphabeta+h_alphabeta[alpha, beta, alpha, beta]
    # ### Second/third terms of eq 26
    # for alpha in range(4):
    #     for i in range(4):
    #         for j in range(4):
    #             R_alpha=R_alpha+h_alpha[alpha, i, j]
    for alpha in range(4):
        for beta in range(4):
            R_alpha=R_alpha+h_alpha[alpha, alpha, beta]
            R_alpha=R_alpha+h_alpha[alpha, beta, alpha]
    ### 4th term of eq 26
    for i in range(4):
        for j in range(4):
            R_og=R_og+h_og[i, j]
    ### 5th
    for alpha in range(4):
        R_2_betabeta=R_2_betabeta+np.trace(h_alphabeta[alpha, alpha])
    ### 6th
    for alpha in range(4):
        R_betabeta=R_betabeta+np.trace(h_alpha[alpha])
    ### 8th
    R_og_betabeta=np.trace(h_og)
    Ricci_scalar=R_alphabeta-R_alpha+R_og-R_2_betabeta+2*R_betabeta-R_og_betabeta
    return Ricci_scalar

### First order approximation of the plaquette. 
### I put the Delta inside the summation, so that I wouldn't have to define another function and restructure gauge_latticeqcd
### Returns only the deltaP, a 4 vector. I still have to contract with the Delta (deformation), which is done in deltaSEH
def Plaq_approx(self, t, x, y, z, matrices, SPrime):
    matrices_length = len(matrices)
    Plaquette_approximations=[0., 0., 0., 0.]
    Action_1=0
    coords=[t, x, y, z]
    for ro in range(4):
        coords[ro]+=1
        for mu in range(4):
            ### This part uses the code from gauge_latticeqcd
            r = np.random.randint(0, matrices_length) 
            matrix = matrices[r]
            staple=self.dS_staple(coords[0], coords[1], coords[2], coords[3], mu)
            link=self.U[coords[0], coords[1], coords[2], coords[3], mu, :, :]
            updated_link = np.dot(matrix, self.U[coords[0], coords[1], coords[2], coords[3], mu, :, :])
            # Plaquette_approximations[ro]=(-1 / 3.0 / self.u0) * np.real(np.trace(np.dot( (updated_link - link), staple)))
            Lqcd_nu=(1 - ((1 / 3.0 / self.u0) * np.real(np.trace(np.dot( (updated_link), staple)))))  ####### Expanded the difference into two lines
            Lqcd_old=(1 - ((1 / 3.0 / self.u0) * np.real(np.trace(np.dot( (link), staple)))))
            Action_1=Action_1+SPrime[t, x, y, z, ro]*Lqcd_nu-self.SP[t, x, y, z, ro]*Lqcd_old
        coords[ro]-=1
    return Action_1