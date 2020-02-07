import matplotlib.pyplot as plt
import numpy as np


'''Setting up reaction scheme
        inital conditions, rate constants, and stoichiometry matrices
'''
#defining initial conditions for species in concentration
#TODO:these need to be convert to number of molecules
C2 = [0.015*(6*10**3)]
CP = [0.015*(6*10**3)]
pM = [0.001*(6*10**3)]
M  = [0.001*(6*10**3)]
Y  = [0.005*(6*10**3)]
YP = [0.005*(6*10**3)]

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
#stoichiometry of total reaction
V_tot =V_r-V_p
# initial concentration
X0 = np.array([Y[0], pM[0], M[0], YP[0], C2[0], CP[0]])

def react(propensity):
    '''Works out which reaction has occurred
    '''
    rand2 = np.random.rand()
    num =0
    for i in propensity:
        if rand2 <=  i:
            num = i 
            break
        
    return propensity.index(num)

# reaction rates
time = [0]
while time[-1] < 1: 
    # rate of each reaction
    k = np.array([k1aaCT, k2, k3CT/(C2[-1] + CP[-1] +pM[-1] + M[-1]), 
    (k4prime+k4*((M[-1]/CT[-1])**2)), k5tilP, k6, k7, k8tilP, k9])
    # current amount of each species evoved in reactions 1 - 9
    species= np.array([1, Y[-1], CP[-1], pM[-1], M[-1], M[-1], YP[-1], C2[-1], CP[-1]])
    # Calculate The total reaction rate
    rates = k*species
    Rtot = rates.sum()
    # Calculate probability of each reaction occurring 
    ProbReact = rates/Rtot

    # distributed propensities
    tot = 0.0
    propensity = []
    for prob in ProbReact:
        tot = tot + prob
        propensity.append(tot)
     # Generate random numbers r1,r2 uniformly distributed in (0,1)
    rand1 = np.random.rand()

    # Compute the time until the next reaction takes place.
    tau = (1.0/ Rtot)*np.log(float(1.0/rand1))
    print('tau: ', tau)
    time.append(time[-1] +tau)
    print('time: ', time[-1])

    #working out which reaction has occured
    reaction = react(propensity)
    #print('reaction choosen: ', reaction)

    # array of current number of molecules
    #molc = np.array([Y[-1], pM[-1], M[-1], YP[-1], C2[-1], CP[-1]])

    #molc = molc - V_r[reaction] + V_p[reaction]

    molc = [Y[-1], pM[-1], M[-1], YP[-1], C2[-1], CP[-1]]
    for i in range(0,len(molc),1):
        molc[i]= molc[i] + V_tot[reaction][i]

    Y.append(molc[0])
    pM.append(molc[1])
    M .append(molc[2])
    YP.append(molc[3])
    C2.append(molc[4])
    CP.append(molc[5])

plt.figure(figsize=(18,10))
plt.title("MPF - SSA")
plt.xlabel("Time (min)")
plt.ylabel("Species Population")
Ynp = np.asarray(Y)
YPnp = np.asarray(YP)
pMnp = np.asarray(pM)
Mnp = np.asarray(M)
CPnp = np.asarray(CP)
C2np = np.asarray(C2)
YT = Ynp + YPnp + pMnp + Mnp
#print('YT: ', YT)
CT = C2np + CPnp + pMnp + Mnp
#print('CT: ', CT)
ytct = YT/CT
myt = Mnp/YT
timenp = np.asarray(time)
#print('myt: ', myt)
#print('ytct: ', ytct)
plt.plot(time ,ytct,'r', label='[YT]/[CT]')
plt.plot(time ,myt,'y', label='[M]/[YT]')

plt.plot([0],[11])
plt.legend(loc='best')
plt.show()

