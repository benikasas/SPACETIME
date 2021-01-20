import numpy as np
from generate import betas, Nt, Nx, Ny, Nz, action, epsilon, u0
import tools_v1
from gauge_latticeqcd import


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
def fluctuations(Nfluc=1000, epsilon=0.2):
    disturbances=[]
    for i in range(Nfluc):
        X=Delta_gen(epsilon)
        disturbances.append(X)
    return disturbances

class spacetime():
    def __init__(self, Nt, Nx, Ny, Nz, beta, u0, U=None):
        if None == U:
            # Creates a spacetime grid, and at every point theres 4 different directions to be perturbed
            U=np.zeros((14, 14, 14, 14, 4))
        # convert to numpy arrays -> significant speed up
        self.U = np.array(U)
        ## Unsure about the necesity of this
        # self.beta = beta
        # self.u0 = u0
        # self.Nx = Nx
        # self.Ny = Ny
        # self.Nz = Nz
        # self.Nt = Nt
    
    ## A function to disturb the initialized empty spacetime array by adding random numbers to U via fluctuations fn.
    ## Decided to put it as a separate function for cleanliness
    def disturbinator(self, disturbances, U):
        U_prime=self.U
        for t in range(Nt):
            for x in range(Nx):
                for y in range(Ny):
                    for z in range(Nz):
                        r=np.random.randint(0, Nfluc)
                        U_prime[t, x, y, z, :]=self.U[t, x, y, z, :] + disturbances[r]
                        ########FIND DELTA S

                        ## Sets the value of our original 
                        # if (np.exp(-1. * dS) > np.random.uniform(0, 1)):
                        #     self.U=U_prime
        U=U_prime
        return self.U


    ## First order approximation in a direction
    ## Make a for loop in the future
    def first_approx_tool(self):
        # directions=[Nt, Nx, Ny, Nz]
        # for i in directions:
        #     U_og=trial.U(Nt, Nx, Ny, Nz)
        #     i=i+1
        #     U_prime_t=trial.U(Nt, Nx, Ny, Nz)       ## Might actually add to the value but might not move forward in direction
        #     partial_approx_t=(U_prime_t-U_og)/lattice_space
        #     approx+=partial_approx_t
        U_og=self.U(Nt, Nx, Ny, Nz)
        U_t=self.U(Nt+1, Nx, Ny, Nz)
        U_x=self.U(Nt, Nx+1, Ny, Nz)
        U_y=self.U(Nt, Nx, Ny+1, Nz)
        U_z=self.U(Nt, Nx, Ny, Nz+1)        
        diff_t=(U_t-U_og)/lattice_space
        diff_x=(U_x-U_og)/lattice_space
        diff_y=(U_y-U_og)/lattice_space
        diff_z=(U_z-U_og)/lattice_space
        total_diff=diff_t+diff_x+diff_y+diff_z
        return total_diff
    
    def Ricci():
    def deltaSEH():

    
def generator():
    U=spacetime(Nt, Nx, Ny, Nz, beta, u0)
    disturbances=fluctuations()
    U_prime=U.disturbinator(disturbances, U)
    # first_approx=first_approx_tool()

value=generator()



