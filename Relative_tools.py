import numpy as np
from generate import betas, Nt, Nx, Ny, Nz, action, epsilon, u0
import tools_v1


Nfluc=1000
beta=betas[0]
lattice_space=tools_v1.fn_a(beta)

## A function to generate a single disturbance in Spacetime
## Need to adjust the magnitude of epsilon
def Delta_gen(epsilon):
    Delta=[0., 0., 0., 0.,]
    for i in range(len(Delta)): 
        Delta[i]=np.random.uniform(-epsilon, epsilon)
    return Delta

## A function to collect a large number of Deltas
def fluctuations(Nfluc,epsilon):
    disturbance=[]
    for i in range(Nfluc):
        X=Delta_gen(epsilon)
        disturbance.append(X)
    return disturbance

class spacetime():
    def __init__(self, Nt, Nx, Ny, Nz, beta, u0, U=None):
        if None == U:
            # Creates a spacetime grid, and at every point theres 4 
            U = [[[[np.zeros(4) for z in range(Nz)] for y in range(Ny)] for x in range(Nx)] for t in range(Nt)]
        # convert to numpy arrays -> significant speed up
        self.U = np.array(U)
        
        ## Unsure about the necesity of this
        # self.beta = beta
        # self.u0 = u0
        # self.Nx = Nx
        # self.Ny = Ny
        # self.Nz = Nz
        # self.Nt = Nt
    
    ## First order approximation in a direction
    def first_approx():
        directions=[Nt, Nx, Ny, Nz]
        for i in directions:
            U_og=trial.U(Nt, Nx, Ny, Nz)
            i=i+1
            U_prime_t=trial.U(Nt, Nx, Ny, Nz)       ## Might actually add to the value but might not move forward in direction
            partial_approx_t=(U_prime_t-U_og)/lattice_space
            approx+=partial_approx_t

value=fluctuations(Nfluc, epsilon)
A=spacetime(Nt, Nx, Ny, Nz, beta, u0)
print(A)