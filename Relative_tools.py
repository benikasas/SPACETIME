import numpy as np
#from generate import Nt, Nx, Ny, Nz, action, epsilon, u0
import tools_v1
import mpmath as mp
# import decimal as dp
# import sympy as sp
### The above modules do not work or at least very tricky to implement for arbitrary precision float calculations    



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


####################################EDITED KAPPA
### Dimensionless Kappa calculator
def kappa_calc(aa):
    kappa=7.6*10**37*aa*aa  ### This was found using the kappa value (from constants) multiplied by lattice spacing^2
    return kappa

## A function to generate a single disturbance in Spacetime
## Need to adjust the magnitude of epsilon
## Still too big??
def Delta_gen(epsilon, magnitude_1):
    Delta=[0., 0., 0., 0.,]
    magnitude_2=epsilon/(magnitude_1)     ### Specify the magnitude of deformations
    for i in range(len(Delta)): 
        Delta[i]=np.random.uniform(-magnitude_2, magnitude_2)
    return Delta

## First order approximations of spacetime deformations in all directions
## Make a for loop in the future
## Could be combined with inv_jack
## This finds the difference between the deformation at (x') and (x'- a * nu) where nu is the direction vector
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
def inv_Jack(SP, t, x, y, z):
    mp.mp.dps = 45
    # jack=mp.matrix(4, 4)
    SP_og=SP[t, x, y, z, :]
    SP_t=SP[t+1, x, y, z, :]
    SP_x=SP[t, x+1, y, z, :]
    SP_y=SP[t, x, y+1, z, :]
    SP_z=SP[t, x, y, z+1, :]
    # diff_t=SP_t-SP_og
    # diff_x=SP_x-SP_og
    # diff_y=SP_y-SP_og
    # diff_z=SP_z-SP_og
    jack_1=np.zeros((4,4), dtype='double')
    jack_1[0, :]=(SP_t-SP_og)
    jack_1[1, :]=(SP_x-SP_og)
    jack_1[2, :]=(SP_y-SP_og)
    jack_1[3, :]=(SP_z-SP_og)
    jack=mp.matrix(jack_1)

    # for i in range(4):
    #     jack[0, i]=mp.mpf(SP_t[i]-SP_og[i])
    # for i in range(4):
    #     jack[1, i]=mp.mpf(SP_x[i]-SP_og[i])
    # for i in range(4):
    #     jack[2, i]=mp.mpf(SP_y[i]-SP_og[i])
    # for i in range(4):
    #     jack[3, i]=mp.mpf(SP_z[i]-SP_og[i])
    for alpha in range(4):
        jack[alpha, alpha]=mp.fadd(1, jack[alpha, alpha], dps=45)

    jack=jack.T
    jack=jack**-1
    # jack=np.transpose(jack)     ### I am pretty sure this is needed, but I might be wrong
    # jack=np.linalg.inv(jack)
    # print(jack)
    return jack

### Finds the h matrix for the S_Prime matrix with given coordinates
### This is found using equation 23
### The only non - zero elements of the sum will be the ones where (nu') = (mu') because of the flat spacetime metric
def h_matrix_producinator(SPrime, coords):
    t=coords[0]
    x=coords[1]
    y=coords[2]
    z=coords[3]
    jack=inv_Jack(SPrime, t, x, y, z)   ### Returns the Jacobian with the given coordinates in the given spacetime matrix
    ### Diagonal terms
    ### I am pretty sure the indices below are correct, though I might be wrong
    ### Need to check it again
    # ### In the case that the indices are not correct, I can undo the transpose in inv_jack
    # h_matrix=np.zeros((4,4))
    # h_matrix[0, 0]=jack[0, 0]*jack[0, 0]+jack[1, 0]*jack[1, 0]+jack[2, 0]*jack[2, 0]+jack[3, 0]*jack[3, 0]
    # h_matrix[1, 1]=jack[0, 1]*jack[0, 1]+jack[1, 1]*jack[1, 1]+jack[2, 1]*jack[2, 1]+jack[3, 1]*jack[3, 1]
    # h_matrix[2, 2]=jack[0, 2]*jack[0, 2]+jack[1, 2]*jack[1, 2]+jack[2, 2]*jack[2, 2]+jack[3, 2]*jack[3, 2]
    # h_matrix[3, 3]=jack[0, 3]*jack[0, 3]+jack[1, 3]*jack[1, 3]+jack[2, 3]*jack[2, 3]+jack[3, 3]*jack[3, 3]

    # ### I am pretty sure the matrix is going to be symmetric.
    # h_matrix[0, 1]=jack[0, 0]*jack[0, 1]+jack[1, 0]*jack[1, 1]+jack[2, 1]*jack[2, 1]+jack[3, 0]*jack[3, 1]
    # h_matrix[0, 2]=jack[0, 0]*jack[0, 2]+jack[1, 0]*jack[1, 2]+jack[2, 2]*jack[2, 2]+jack[3, 0]*jack[3, 2]
    # h_matrix[0, 3]=jack[0, 0]*jack[0, 3]+jack[1, 0]*jack[1, 3]+jack[2, 3]*jack[2, 3]+jack[3, 0]*jack[3, 3]
    # h_matrix[1, 0]=h_matrix[0, 1]
    # h_matrix[2, 0]=h_matrix[0, 2]
    # h_matrix[3, 0]=h_matrix[0, 3]

    # h_matrix[1, 2]=jack[0, 1]*jack[0, 2]+jack[1, 1]*jack[1, 2]+jack[2, 1]*jack[2, 2]+jack[3, 1]*jack[3, 2]
    # h_matrix[1, 3]=jack[0, 1]*jack[0, 3]+jack[1, 1]*jack[1, 3]+jack[2, 1]*jack[2, 3]+jack[3, 1]*jack[3, 3]
    # h_matrix[2, 1]=h_matrix[1, 2]
    # h_matrix[3, 1]=h_matrix[1, 3]

    # h_matrix[2, 3]=jack[0, 2]*jack[0, 3]+jack[1, 2]*jack[1, 3]+jack[2, 2]*jack[2, 3]+jack[3, 2]*jack[3, 3]
    # h_matrix[3, 2]=h_matrix[2, 3]
    # storage=h_matrix


    h_matrix=mp.matrix(4, 4)
    h_matrix[0, 0]=mp.fadd(mp.fmul(jack[0, 0], jack[0, 0]), mp.fadd(mp.fmul(jack[1, 0], jack[1, 0]), mp.fadd(mp.fmul(jack[2, 0], jack[2, 0]), mp.fmul(jack[3, 0], jack[3, 0]))))
    h_matrix[1, 1]=mp.fadd(mp.fmul(jack[0, 1], jack[0, 1]), mp.fadd(mp.fmul(jack[1, 1], jack[1, 1]), mp.fadd(mp.fmul(jack[2, 1], jack[2, 1]), mp.fmul(jack[3, 1], jack[3, 1]))))
    h_matrix[2, 2]=mp.fadd(mp.fmul(jack[0, 2], jack[0, 2]), mp.fadd(mp.fmul(jack[1, 2], jack[1, 2]), mp.fadd(mp.fmul(jack[2, 2], jack[2, 2]), mp.fmul(jack[3, 2], jack[3, 2]))))
    h_matrix[3, 3]=mp.fadd(mp.fmul(jack[0, 3], jack[0, 3]), mp.fadd(mp.fmul(jack[1, 3], jack[1, 3]), mp.fadd(mp.fmul(jack[2, 3], jack[2, 3]), mp.fmul(jack[3, 3], jack[3, 3]))))
    h_matrix[0, 0]=mp.fadd(-1, h_matrix[0, 0], dps=50)       ################ADDING 1 because its in the metric??
    h_matrix[1, 1]=mp.fadd(-1, h_matrix[1, 1], dps=50)
    h_matrix[2, 2]=mp.fadd(-1, h_matrix[2, 2], dps=50)
    h_matrix[3, 3]=mp.fadd(-1, h_matrix[3, 3], dps=50)

    h_matrix[0, 1]=mp.fadd(mp.fmul(jack[0, 0], jack[0, 1]), mp.fadd(mp.fmul(jack[1, 0], jack[1, 1]), mp.fadd(mp.fmul(jack[2, 0], jack[2, 1]), mp.fmul(jack[3, 0], jack[3, 1]))))
    h_matrix[0, 2]=mp.fadd(mp.fmul(jack[0, 0], jack[0, 2]), mp.fadd(mp.fmul(jack[1, 0], jack[1, 2]), mp.fadd(mp.fmul(jack[2, 0], jack[2, 2]), mp.fmul(jack[3, 0], jack[3, 2]))))
    h_matrix[0, 3]=mp.fadd(mp.fmul(jack[0, 0], jack[0, 3]), mp.fadd(mp.fmul(jack[1, 0], jack[1, 3]), mp.fadd(mp.fmul(jack[2, 0], jack[2, 3]), mp.fmul(jack[3, 0], jack[3, 3]))))
    h_matrix[1, 0]=h_matrix[0, 1]
    h_matrix[2, 0]=h_matrix[0, 2]
    h_matrix[3, 0]=h_matrix[0, 3]
    h_matrix[1, 2]=mp.fadd(mp.fmul(jack[0, 1], jack[0, 2]), mp.fadd(mp.fmul(jack[1, 1], jack[1, 2]), mp.fadd(mp.fmul(jack[2, 1], jack[2, 2]), mp.fmul(jack[3, 1], jack[3, 2]))))
    h_matrix[1, 3]=mp.fadd(mp.fmul(jack[0, 1], jack[0, 3]), mp.fadd(mp.fmul(jack[1, 1], jack[1, 3]), mp.fadd(mp.fmul(jack[2, 1], jack[2, 3]), mp.fmul(jack[3, 1], jack[3, 3]))))
    h_matrix[2, 1]=h_matrix[1, 2]
    h_matrix[3, 1]=h_matrix[1, 3]
    h_matrix[2, 3]=mp.fadd(mp.fmul(jack[0, 2], jack[0, 3]), mp.fadd(mp.fmul(jack[1, 2], jack[1, 3]), mp.fadd(mp.fmul(jack[2, 2], jack[2, 3]), mp.fmul(jack[3, 2], jack[3, 3]))))
    h_matrix[3, 2]=h_matrix[2, 3]

    A=np.zeros((4, 4))
    for i in range(4):
        for j in range(4):
            A[i, j]=float(h_matrix[i, j])
    ### Converts h_matrix into a numpy array, faster, easier, 
    # h_matrix[0, 0]=jack[1][0]*jack[1][0]+jack[2][0]*jack[2][0]+jack[3][0]*jack[3][0]
    # h_matrix[1, 1]=jack[0][1]*jack[0][1]+jack[2][1]*jack[2][1]+jack[3][1]*jack[3][1]
    # h_matrix[2, 2]=jack[0][2]*jack[0][2]+jack[1][2]*jack[1][2]+jack[3][2]*jack[3][2]
    # h_matrix[3, 3]=jack[0][3]*jack[0][3]+jack[1][3]*jack[1][3]+jack[2][3]*jack[2][3]

    return A

### Collects surrounding h matrices for given coordinates into a 4x4 matrix 
### Returns the h matrices separated by two (not necessarily different) unit vectors
### Returns also the h matrices separated by 1 unit vector from the origin
### Lastly, returns the h matrix for the original coordinates
### Could split into 3 different functions, but this saves 1 line of code I think
def h_matrix_collectinator(SPrime, t, x, y, z):
    h_alphabeta=[[np.identity(4) for beth in range(4)] for alph in range(4)]    ### For h matrices where we add two lattice spacings from the origin
    h_alpha=[np.identity(4) for alph in range(4)]   ### For h matrices separated from origin by single lattice spacing
    coords=[t, x, y, z]
    for alpha in range(4):      ### Used for first unit vector
        coords[alpha]+=1    ### Shifts one coordinate by 1
        for beta in range(4):   ### Used for the second unit vector
            coords[beta]+=1     ### Shifts another coordinate by 1
            h_alphabeta[alpha][beta]=h_matrix_producinator(SPrime, coords)
            coords[beta]-=1     ### Returns the second coordinate into the original position
        h_alpha[alpha]=h_matrix_producinator(SPrime, coords)
        coords[alpha]-=1        ### Returns the first coordinate
    h_og=h_matrix_producinator(SPrime, coords)      ### The h matrix for origin
    return h_alphabeta, h_alpha, h_og


### Function to find the Ricci scalar curvature at a given point as given in eq. 26
def Ricci(SPrime, t, x, y, z):
    h_alphabeta, h_alpha, h_og = h_matrix_collectinator(SPrime, t, x, y, z)
    h_alphabeta=np.array(h_alphabeta)
    h_alpha=np.array(h_alpha)
    h_og=np.array(h_og)
    R_alphabeta=0.       ### Used for first term of the sum (1st partial sum)
    R_alpha=0.           ### 2nd term and 3rd term summed up
    R_og=0.              ### 4th term
    R_2_betabeta=0.      ### 5th term
    R_betabeta=0.        ### 6th term
    R_og_betabeta=0.     ### 7th term
    Ricci_scalar=0.      ### The final sum
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
            R_alpha=R_alpha+h_alpha[alpha, beta, alpha]     ### I cut some corners, h_(alpha, beta)(x+a*beta)=h_(beta, alpha)(x+a*beta) as we sum over all alphas and betas
    ### 4th term of eq 26
    for i in range(4):
        for j in range(4):
            R_og=R_og+h_og[i, j]
    ### 5th
    for alpha in range(4):
        R_2_betabeta=R_2_betabeta+np.trace(h_alphabeta[alpha, alpha])   ### For every alpha, the h_(beta, beta) will return the trace when summed over beta
    ### 6th
    for alpha in range(4):
        R_betabeta=R_betabeta+np.trace(h_alpha[alpha])      ### Same as above
    ### 8th
    R_og_betabeta=np.trace(h_og)
    Ricci_scalar=R_alphabeta-R_alpha+R_og-R_2_betabeta+2*R_betabeta-R_og_betabeta   ### Finds the total sum
    return Ricci_scalar

### First order approximation of the plaquette. 
### I put the Delta inside the summation of eq. 22, so that I wouldn't have to define another function and restructure gauge_latticeqcd
def Plaq_approx(self, t, x, y, z, matrices, SPrime):
    matrices_length = len(matrices)
    Plaquette_approximations=[0., 0., 0., 0.]
    Action_2=0
    coords=[t, x, y, z]
    for ro in range(4):
        coords[ro]+=1       ### Pushes the coordinates by one lattice link to each direction
        for mu in range(4): ### Maybe its this mu that causes such a significant difference??
            ### This part uses the code from gauge_latticeqcd
            r = np.random.randint(0, matrices_length) 
            matrix = matrices[r]
            staple=self.dS_staple(coords[0], coords[1], coords[2], coords[3], mu)
            link=self.U[coords[0], coords[1], coords[2], coords[3], mu, :, :]
            updated_link = np.dot(matrix, self.U[coords[0], coords[1], coords[2], coords[3], mu, :, :])
            # Plaquette_approximations[ro]=(-1 / 3.0 / self.u0) * np.real(np.trace(np.dot( (updated_link - link), staple)))
            Lqcd_nu=(1 - ((1 / 3.0 / self.u0) * np.real(np.trace(np.dot( (updated_link), staple)))))  ####### Expanded the difference into two lines
            Lqcd_old=(1 - ((1 / 3.0 / self.u0) * np.real(np.trace(np.dot( (link), staple)))))
            Action_2=Action_2+SPrime[t, x, y, z, ro]*Lqcd_nu-self.SP[t, x, y, z, ro]*Lqcd_old
        coords[ro]-=1
    return Action_2

