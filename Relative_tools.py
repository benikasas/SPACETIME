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

## First order approximation in a direction
## Make a for loop in the future
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
        total_diff=diff_t+diff_x+diff_y+diff_z
        total_sum=0
        for i in range(len(total_diff)):
            total_sum=total_sum+total_diff[i]
            print(i, ' term: ', total_diff[i])
            print('Total sum after ', i, ' ', total_sum)
        return 1

### First order approximation of the Jacobian inverse
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
    jack=np.linalg.inv(jack)
    return jack