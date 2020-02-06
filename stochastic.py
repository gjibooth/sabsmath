import gillespy2
import sys, os
import matplotlib.pyplot as plt
import numpy


class MPF(gillespy2.Model):
       def __init__(self, parameter_values=None):
              # initialize Model
              gillespy2.Model.__init__(self, name="MPF")

              ''' 
              PARAMETERS:
              These are constant values relevant to the system, such as reaction kinetic rates.
  
              name: user defined name for reference
              expression: constant value
              '''
              k1aaCT = gillespy2.Parameter(name='k1aaCT', expression=0.015)
              k2 = gillespy2.Parameter(name='k2', expression=0.0)
              k3CT = gillespy2.Parameter(name='k3CT', expression=200.0)
              k4 = gillespy2.Parameter(name='k4', expression=180.0)
              k4prime = gillespy2.Parameter(name='k4prime', expression=0.018)
              k5tilP = gillespy2.Parameter(name='k5tilP', expression=0.0)
              k6 = gillespy2.Parameter(name='k6', expression=1.0)
              k7 = gillespy2.Parameter(name='k7', expression=0.6)
              k8tilP = gillespy2.Parameter(name='k8tilP', expression=50.0)
              k9 = gillespy2.Parameter(name='k9', expression=10.0)


              # Add parameters to the model
              self.add_parameter([k1aaCT, k2, k3CT, k4, k4prime, k5tilP, k6, k7, k8tilP, k9])

              '''
              SPECIES:
              These can be anything that participates in or is produced by a reaction channel.
  
              name: A user defined name for the species
              initial_value: value/population count of species at start of simulation
              '''
              #cdc2
              C2 = gillespy2.Species(name='C2', initial_value=301)
              #cdc2p
              CP = gillespy2.Species(name='CP', initial_value=120)
              #cdc2p inactive
              pM = gillespy2.Species(name='pM', initial_value=0)
              #active dimer cyc-p-cdc2
              M = gillespy2.Species(name='D', initial_value=0)
              #cyc
              Y = gillespy2.Species(name='Y', initial_value=120)
              #cyc-p
              YP = gillespy2.Species(name='YP', initial_value=0)
              #aminoacids
              aa = gillespy2.Species(name= 'aa', initial_value=1)



              # Add species to the model
              self.add_species([C2, CP, pM, M, Y, YP, aa])

              '''Reactions:
              These are the reaction channels causing the system to change over time
  
              name: a user defined name for the reaction
              reactants: dictionary with participant reactants as keys, and consumed per reaction as value.
              products: dictionary with reaction products as keys, and number formed per reaction as value.
              rate: parameter rate constant to be applied to the propensity of this reaction firing
              propensity_function: can be used instead of rate in order to declare a custom propensity function in string format'''
              r1 = gillespy2.Reaction(name="r1", reactants={aa: 1}, products={Y: 1},
                                      rate=k1aaCT)

              r2 = gillespy2.Reaction(name="r2", reactants={Y: 1}, products={aa: 1},
                                      rate=k2)

              r3 = gillespy2.Reaction(name="r3", reactants={Y: 1}, products={pM: 1},
                                      rate=k3CT)

              r4 = gillespy2.Reaction(name="r4", reactants={pM: 1}, products={M: 1},
                                      rate=k4)

              r5 = gillespy2.Reaction(name="r5", reactants={M: 1}, products={pM: 1},
                                      rate=k5tilP)

              r6 = gillespy2.Reaction(name="r6", reactants={M: 1}, products={C2: 1, YP: 1},
                                      rate=k6)

              r7 = gillespy2.Reaction(name="r7", reactants={YP: 1}, products={aa: 1},
                                      rate=k7)

              r8 = gillespy2.Reaction(name="r8", reactants={C2: 1}, products={CP: 1},
                                      rate=k8tilP)

              r9 = gillespy2.Reaction(name="r9", reactants={CP: 1}, products={C2: 1},
                                      rate=k9)

              # Add reactions to the model
              self.add_reaction([r1, r2, r3, r4, r5, r6, r7, r8, r9])

              # Set timespan of model
              self.timespan(numpy.linspace(0, 100, 101))


# Instantiate your model
model = MPF()

'''Run a stochastic simulation on the model and store results to a variable.  
If a solver is not explicitly declared, GillesPy2 will select the direct SSA method'''
results = model.run()

results.plot()
plt.show()
print('here')