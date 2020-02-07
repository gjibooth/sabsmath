import matplotlib.pyplot as plt
import numpy as np


'''Setting up reaction scheme
        inital conditions, rate constants, and stoichiometry matrices
'''
#defining initial conditions for species in concentration
#TODO:these need to be convert to number of molecules
C2 = [0.015]
CP = [0.015]
pM = [0.001]
M  = [0.001]
Y  = [0.005]
YP = [0.005]
t = [0]
tf = [30]

YT = [Y[0] + YP[0] + pM[0] + M[0]]
CT = [C2[0] + CP[0] + pM[0] + M[0]]

# Rate Constants
k1aaCT  = 0.015
k2      = 0
k3CT    = 200
k4      = 180 # adjustable
k4prime = 0.018
k5tilP  = 0
k6      = 1 # adjustable
k7      = 0.6
k8tilP  = 50 # >> k9
k9      = 10 # >> k6

# reactant stoichiometry
#                Y,pM, M,YP,C2,CP 
V_r = np.array([[0, 0, 0, 0, 0, 0], 
                [1, 0, 0, 0, 0, 0], 
                [0, 0, 0, 0, 0, 1],
                [0, 1, 0, 0, 0, 0],
                [0, 0, 1, 0, 0, 0],
                [0, 0, 1, 0, 0, 0],
                [0, 0, 0, 1, 0, 0],
                [0, 0, 0, 0, 1, 0],
                [0, 0, 0, 0, 0, 1]])
# product stoichiometry
#                Y,pM, M,YP,C2,CP 
V_p = np.array([[1, 0, 0, 0, 0, 0], 
                [0, 0, 0, 0, 0, 0], 
                [0, 1, 0, 0, 0, 0],
                [0, 0, 1, 0, 0, 0],
                [0, 1, 0, 0, 0, 0], 
                [0, 0, 0, 1, 1, 0],
                [0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 1],
                [0, 0, 0, 0, 1, 0]])
# initial concentration
X0 = np.array([Y[0], pM[0], M[0], YP[0], C2[0], CP[0]])

# number steps to be used in 
timespan = np.linspace(0, 100, 101)
# reaction rates
time = 0
for t in timespan:
    # rate of each reaction
    k = np.array([k1aaCT, k2, k3CT/(C2[-1] + CP[-1] +pM[-1] + M[-1]), 
    (k4prime+k4*((M[-1]/CT[-1])**2)), k5tilP, k6, k7, k8tilP, k9])
    # current amount of each species evoved in reactions 1 - 9
    species= np.array([1, Y[-1], CP[-1], pM[-1], M[-1], M[-1], YP[-1], C2[-1], CP[-1]])
    # Calculate The total reaction rate
    rates = k*species
    Rtot= rates.sum()
    # Calculate probability of each reaction occurring 
    ProbReact = rates/Rtot
    
    #TODO: Calculate weighted probability of each reaction occurring i.e propensity functions
    
    # Generate random numbers r1,r2 uniformly distributed in (0,1)
    r1 = np.random.rand()
    r2 = np.random.rand()

# if using tau leaping
#TODO:calculate alpha
alpha = 1.0

# Compute the time until the next reaction takes place.
tau = (1.0/alpha)*np.log(float(1.0/r1))
time = time +tau