import numpy as np

Nfluc=1000
epsilon=0.2


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

for i in range(10):
    print(i)
    i=i+1
    print(i)

value=fluctuations(Nfluc, epsilon)