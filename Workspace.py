import numpy as np
from gauge_latticeqcd import matrix_su2, matrix_su3, create_su3_set, lattice

Nt=5
Nz=5
Ny=5
Nx=5
U = [[[[[np.identity(3, dtype='complex128') for mu in range(4)] for z in range(Nz)] for y in range(Ny)] for x in range(Nx)] for t in range(Nt)]
print(U[1][1][1][1])