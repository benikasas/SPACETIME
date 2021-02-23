import numpy as np
# from gauge_latticeqcd import matrix_su2, matrix_su3, create_su3_set, lattice
# from generate import *
# import tools_v1 as tool
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
    h_matrix=np.zeros((4,4))
    h_matrix[0, 0]=jack[0, 0]*jack[0, 0]+jack[1, 0]*jack[1, 0]+jack[2, 0]*jack[2, 0]+jack[3, 0]*jack[3, 0]
    h_matrix[1, 1]=jack[0, 1]*jack[0, 1]+jack[1, 1]*jack[1, 1]+jack[2, 1]*jack[2, 1]+jack[3, 1]*jack[3, 1]
    h_matrix[2, 2]=jack[0, 2]*jack[0, 2]+jack[1, 2]*jack[1, 2]+jack[2, 2]*jack[2, 2]+jack[3, 2]*jack[3, 2]
    h_matrix[3, 3]=jack[0, 3]*jack[0, 3]+jack[1, 3]*jack[1, 3]+jack[2, 3]*jack[2, 3]+jack[3, 3]*jack[3, 3]

    h_matrix=mp.matrix((4, 4))
    h_matrix[0, 0]=mp.fmul(jack[0, 0], jack[0, 0])+mp.fmul(jack[1, 0], jack[1, 0])+mp.fmul(jack[2, 0], jack[2, 0])+mp.fmul(jack[3, 0], jack[3, 0])
    h_matrix[1, 1]=mp.fmul(jack[0, 1], jack[0, 1])+mp.fmul(jack[1, 1], jack[1, 1])+mp.fmul(jack[2, 1], jack[2, 1])+mp.fmul(jack[3, 1], jack[3, 1])
    h_matrix[2, 2]=mp.fmul(jack[0, 2], jack[0, 2])+mp.fmul(jack[1, 2], jack[1, 2])+mp.fmul(jack[2, 2], jack[2, 2])+mp.fmul(jack[3, 2], jack[3, 2])
    h_matrix[3, 3]=mp.fmul(jack[0, 3], jack[0, 3])+mp.fmul(jack[1, 3], jack[1, 3])+mp.fmul(jack[2, 3], jack[2, 3])+mp.fmul(jack[3, 3], jack[3, 3])
    h_matrix[0, 0]=mp.fadd(-1, h_matrix[0, 0], dps=50)
    h_matrix[1, 1]=mp.fadd(-1, h_matrix[1, 1], dps=50)
    h_matrix[2, 2]=mp.fadd(-1, h_matrix[2, 2], dps=50)
    h_matrix[3, 3]=mp.fadd(-1, h_matrix[3, 3], dps=50)

    # h_matrix[0, 0]=jack[1][0]*jack[1][0]+jack[2][0]*jack[2][0]+jack[3][0]*jack[3][0]
    # h_matrix[1, 1]=jack[0][1]*jack[0][1]+jack[2][1]*jack[2][1]+jack[3][1]*jack[3][1]
    # h_matrix[2, 2]=jack[0][2]*jack[0][2]+jack[1][2]*jack[1][2]+jack[3][2]*jack[3][2]
    # h_matrix[3, 3]=jack[0][3]*jack[0][3]+jack[1][3]*jack[1][3]+jack[2][3]*jack[2][3]


    ### I am pretty sure the matrix is going to be symmetric.
    h_matrix[0, 1]=jack[0, 0]*jack[0, 1]+jack[1, 0]*jack[1, 1]+jack[2, 1]*jack[2, 1]+jack[3, 0]*jack[3, 1]
    h_matrix[0, 2]=jack[0, 0]*jack[0, 2]+jack[1, 0]*jack[1, 2]+jack[2, 2]*jack[2, 2]+jack[3, 0]*jack[3, 2]
    h_matrix[0, 3]=jack[0, 0]*jack[0, 3]+jack[1, 0]*jack[1, 3]+jack[2, 3]*jack[2, 3]+jack[3, 0]*jack[3, 3]
    h_matrix[1, 0]=h_matrix[0, 1]
    h_matrix[2, 0]=h_matrix[0, 2]
    h_matrix[3, 0]=h_matrix[0, 3]


    h_matrix[1, 2]=jack[0, 1]*jack[0, 2]+jack[1, 1]*jack[1, 2]+jack[2, 1]*jack[2, 2]+jack[3, 1]*jack[3, 2]
    h_matrix[1, 3]=jack[0, 1]*jack[0, 3]+jack[1, 1]*jack[1, 3]+jack[2, 1]*jack[2, 3]+jack[3, 1]*jack[3, 3]
    h_matrix[2, 1]=h_matrix[1, 2]
    h_matrix[3, 1]=h_matrix[1, 3]

    h_matrix[2, 3]=jack[0, 2]*jack[0, 3]+jack[1, 2]*jack[1, 3]+jack[2, 2]*jack[2, 3]+jack[3, 2]*jack[3, 3]
    h_matrix[3, 2]=h_matrix[2, 3]
    return h_matrix
