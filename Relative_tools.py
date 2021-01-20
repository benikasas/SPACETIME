import numpy as np
from generate import betas, Nt, Nx, Ny, Nz, action, epsilon, u0
import tools_v1


Nfluc=1000
beta=betas[0]
lattice_space=tools_v1.fn_a(beta)



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

## A function to generate a single disturbance in Spacetime
## Need to adjust the magnitude of epsilon
def Delta_gen(epsilon):
    Delta=[0., 0., 0., 0.,]
    for i in range(len(Delta)): 
        Delta[i]=np.random.uniform(-epsilon, epsilon)
    return Delta

## A function to collect a large number of Deltas
# def fluctuations(Nfluc=1000, epsilon=0.2):
#     disturbances=[]
#     for i in range(Nfluc):
#         X=Delta_gen(epsilon)
#         disturbances.append(X)
#     return disturbances


## First order approximation in a direction
## Make a for loop in the future
def first_approx_tool(SP, t, x, y, z):
        #     i=i+1
        #     U_prime_t=trial.U(Nt, Nx, Ny, Nz)       ## Might actually add to the value but might not move forward in direction
        #     partial_approx_t=(U_prime_t-U_og)/lattice_space
        #     approx+=partial_approx_t
        # for i in range(4):
        #     SP_og=SP[t, x, y, z, :]
        SP_og=SP[t, x, y, z, :]
        SP_t=SP[t+1, x, y, z, :]
        SP_x=SP[t, x+1, y, z, :]
        SP_y=SP[t, x, y+1, z, :]
        SP_z=SP[t, x, y, z+1, :]        
        # diff_t=(SP_t-SP_og)/lattice_space
        # diff_x=(SP_x-SP_og)/lattice_space
        # diff_y=(SP_y-SP_og)/lattice_space
        # diff_z=(SP_z-SP_og)/lattice_space
        ### AAAAAAA
        diff_t=SP_t-SP_og
        diff_x=SP_x-SP_og
        diff_y=SP_y-SP_og
        diff_z=SP_z-SP_og     
        total_diff=diff_t+diff_x+diff_y+diff_z
        # length=np.linalg.norm(total_diff)#### WHAT DO I DOOOOO
        return 1

### First order approximation of the Jacobian inverse
def inv_Jack(SP, t, x, y, z):
    jack=np.zeros((4, 4))
    SP_og=SP[t, x, y, z, :]
    SP_t=SP[t+1, x, y, z, :]
    SP_x=SP[t, x+1, y, z, :]
    SP_y=SP[t, x, y+1, z, :]
    SP_z=SP[t, x, y, z+1, :]
    jack[0][:]=(SP_t-SP_og)
    jack[1][:]=(SP_x-SP_og)
    jack[2][:]=(SP_y-SP_og)
    jack[3][:]=(SP_z-SP_og)
    print(jack)
    return jack



# class spacetime():
#     def __init__(self, Nt, Nx, Ny, Nz, beta, u0, SP=None):
#         if None == SP:
#             # Creates a spacetime grid, and at every point theres 4 different directions to be perturbed
#             SP=np.zeros((14, 14, 14, 14, 4))
#         # convert to numpy arrays -> significant speed up
#         self.SP = np.array(SP)
#         ## Unsure about the necesity of this
#         # self.beta = beta
#         # self.u0 = u0
#         # self.Nx = Nx
#         # self.Ny = Ny
#         # self.Nz = Nz
#         # self.Nt = Nt
    
#     ## A function to disturb the initialized empty spacetime array by adding random numbers to U via fluctuations fn.
#     ## Decided to put it as a separate function for cleanliness
#     def disturbinator(self, disturbances, U):
#         SP_prime=self.SP
#         for t in range(Nt):
#             for x in range(Nx):
#                 for y in range(Ny):
#                     for z in range(Nz):
#                         r=np.random.randint(0, Nfluc)
#                         U_prime[t, x, y, z, :]=self.U[t, x, y, z, :] + disturbances[r]
#                         ########FIND DELTA S

#                         ## Sets the value of our original 
#                         # if (np.exp(-1. * dS) > np.random.uniform(0, 1)):
#                         #     self.U=U_prime
#         SP=SP_prime
#         return self.SP


#     ## First order approximation in a direction
#     ## Make a for loop in the future
#     def first_approx_tool(self):
#         # directions=[Nt, Nx, Ny, Nz]
#         # for i in directions:
#         #     U_og=trial.U(Nt, Nx, Ny, Nz)
#         #     i=i+1
#         #     U_prime_t=trial.U(Nt, Nx, Ny, Nz)       ## Might actually add to the value but might not move forward in direction
#         #     partial_approx_t=(U_prime_t-U_og)/lattice_space
#         #     approx+=partial_approx_t
#         SP_og=self.SP(Nt, Nx, Ny, Nz)
#         SP_t=self.SP(Nt+1, Nx, Ny, Nz)
#         SP_x=self.SP(Nt, Nx+1, Ny, Nz)
#         SP_y=self.SP(Nt, Nx, Ny+1, Nz)
#         SP_z=self.SP(Nt, Nx, Ny, Nz+1)        
#         diff_t=(SP_t-SP_og)/lattice_space
#         diff_x=(SP_x-SP_og)/lattice_space
#         diff_y=(SP_y-SP_og)/lattice_space
#         diff_z=(SP_z-SP_og)/lattice_space
#         total_diff=diff_t+diff_x+diff_y+diff_z
#         return total_diff
    
# def Ricci():


# def deltaSEH():

    
# def generator():
#     U=spacetime(Nt, Nx, Ny, Nz, beta, u0)
#     disturbances=fluctuations()
#     U_prime=U.disturbinator(disturbances, U)
#     # first_approx=first_approx_tool()

# value=generator()



