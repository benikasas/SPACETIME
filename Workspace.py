import numpy as np
# from gauge_latticeqcd import matrix_su2, matrix_su3, create_su3_set, lattice
# from generate import *
# import tools_v1 as tool

# beta=betas[0]
Nt=5
Nz=7
Ny=9
Nx=10
u0=1
# lattice_space
# class spacetime():
#     def __init__(self, Nt, Nx, Ny, Nz, beta, u0, U=None):
#         if None == U:
#             # Creates a spacetime grid, and at every point theres 4 
#             U = [[[[np.zeros(4) for z in range(Nz)] for y in range(Ny)] for x in range(Nx)] for t in range(Nt)]
#         # convert to numpy arrays -> significant speed up
#         self.U = np.array(U)

#         ## Unsure if needed
#         # self.beta = beta
#         # self.u0 = u0
#         # self.Nx = Nx
#         # self.Ny = Ny
#         # self.Nz = Nz
#         # self.Nt = Nt

# # def first_approx():
#     U_og=trial.U(Nt, Nx, Ny, Nz)-trial.U(Nt, Nx, Ny, Nz)
#     U_prime=trial.U(Nt, Nx, Ny, Nz)
#     approx=(U_prime-U_og)/lattice
# Trial=spacetime(Nt, Nx, Ny, Nz, beta, u0)
# print(Trial.U[1, 1, 1, 1])
# partial_sum=0
# jack=np.zeros((4,4))
# jack[1,1]=1
# for i in range(4):
#     for j in range(4):
#         print(jack[i][j])

# SP=np.zeros((14, 14, 14, 14, 4))
# print(len(SP[])

# def inv_Jack(SPrime, t, x, y, z):
#     SP_og=SPrime[t, x, y, z, :]
#     SP_t=SPrime[t+1, x, y, z, :]
#     SP_x=SPrime[t, x+1, y, z, :]
#     SP_y=SPrime[t, x, y+1, z, :]
#     SP_z=SPrime[t, x, y, z+1, :]
#     jack=np.identity(4)
#     jack[0][:]=np.add(jack[0][:], (SP_t-SP_og))
#     jack[1][:]=np.add(jack[1][:], (SP_x-SP_og))
#     jack[2][:]=np.add(jack[2][:], (SP_y-SP_og))
#     jack[3][:]=np.add(jack[3][:], (SP_z-SP_og))
#     jack=np.linalg.inv(jack)
#     return jack

# def Ricci(t, x, y, z):
#     array=[t, x, y, z]
#     for alpha in range(4):
#         for beta in range(4):
#             array[alpha]=array[alpha]+1
#             array[beta]=array[beta]+1
#             print(array)
#             array[alpha]-=1
#             array[beta]-=1
#             print('space1')
#             if alpha == beta:
#                 print('ayoo')
#         array[alpha]+=1
#         print(array)
#         array[alpha]-=1
#         print('space2')
#     return 1

# Ricci(14, 14, 14, 14)

# def kappa_calca(aa):
#     kappa=7.6*10**37*aa*aa
#     return kappa

# def fn_a(beta):
#     return 0.5 * np.exp(-1.6804 - 1.7331 * (beta - 6.) + 0.7849 * (beta - 6.)**2 - 0.4428 * (beta - 6.)**3)

# aa=fn_a(6.92)
# print(aa)
# print(kappa_calca(aa))