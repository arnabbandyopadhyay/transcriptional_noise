# -*- coding: utf-8 -*-
"""
Created on Sun Nov 21 20:24:59 2021

@author: Arnab
"""

# Numpy is used to set the timespan of the Model.
import numpy as np
import matplotlib.pyplot as plt




# Import the types that'll be needed to define your Model.
from gillespy2.core import (
    Model,
    Species,
    Reaction,
    Parameter
)

from gillespy2.solvers.numpy import (
 	NumPySSASolver,
 	ODESolver,
 	TauLeapingSolver,
 	TauHybridSolver
)


from gillespy2.solvers.cpp import (
 	SSACSolver,
 	ODECSolver,
 	TauLeapingCSolver
)

"""
Your model is declared and configured as a Python class. As such, the name 
between `class` and `(Model):` can be of your choosing.

For this example we'll be modeling a Michaelis-Menten enzymatic reaction,
so lets set the name accordingly.
"""

class MichaelisMenten(Model):
     def __init__(self, pinpd, pinpk, pactpd, pactpk, ppr, pg3p, pmd, pmk, pmr, pd, pk, pr, pgly, parameter_values=None):

            # Intialize the Model with a name of your choosing.
            Model.__init__(self, name="Michaelis_Menten")
            
            """
            Parameters are constant values relevant to the system, such as reaction kinetic rates.
            
            - name: A user defined name for reference.
            - expression: Some constant value.
            """
           
            rate1 = Parameter(name="konD", expression=0.001)
            rate2 = Parameter(name="konK", expression=0.001)
            rate3 = Parameter(name="koffD", expression=1)
            rate4 = Parameter(name="koffK", expression=1)
            rate5 = Parameter(name="ktrsc", expression=0.5)
            rate6 = Parameter(name="ktrsc1", expression=1)
            rate7 = Parameter(name="ktrsc2", expression=1/10)
            rate8 = Parameter(name="ktrsc3", expression=1/10)
            rate9 = Parameter(name="ktrsc4", expression=1/15)
            rate10 = Parameter(name="ktrsl", expression=0.5)
            rate11 = Parameter(name="ktrsl1", expression=0.5)
            rate12 = Parameter(name="ktrsl2", expression=0.5)
            rate13 = Parameter(name="ktrsp_g3p", expression=0)
            rate14 = Parameter(name="ktrsp_gly", expression=0)
            rate15 = Parameter(name="kpf", expression=0.005)
            rate16 = Parameter(name="knf", expression=0.005)
            rate17 = Parameter(name="dmD", expression=0.069)
            rate18 = Parameter(name="dmK", expression=0.069)
            rate19 = Parameter(name="dmR", expression=0.069)
            rate20 = Parameter(name="dpD", expression=0.023)
            rate21 = Parameter(name="dpK", expression=0.023)
            rate22 = Parameter(name="dpR", expression=0.023)
            
            
            # Add the Parameters to the Model.
            self.add_parameter([rate1, rate2, rate3, rate4, rate5, rate6,
                                rate7, rate8, rate9, rate10, rate11, rate12,
                                rate13, rate14, rate15, rate16, rate17, rate18,
                                rate19, rate20, rate21, rate22])
            
            """
            Species can be anything that participates in or is produced by a reaction channel.
            
            - name: A user defined name for the species.
            - initial_value: A value/population count of species at start of simulation.
            """

            inPD_R = Species(name="inPD_R", initial_value=pinpd)
            inPK_R = Species(name="inPK_R", initial_value=pinpk)
            actPD = Species(name="actPD", initial_value=pactpd)
            actPK = Species(name="actPK", initial_value=pactpk)
            PR = Species(name="PR", initial_value=ppr)
            G3P = Species(name="G3P", initial_value=pg3p)
            mD = Species(name="mD", initial_value=pmd)
            mK = Species(name="mK", initial_value=pmk)
            mR = Species(name="mR", initial_value=pmr)
            D = Species(name="D", initial_value=pd)
            K = Species(name="K", initial_value=pk)
            R = Species(name="R", initial_value=pr)
            Gly = Species(name="Gly", initial_value=pgly)
            
            # Add the Species to the Model.
            self.add_species([inPD_R, inPK_R, actPD, actPK, PR, G3P, mD, mK, mR, D, K, R, Gly])
            
            """
            Reactions are the reaction channels which cause the system to change over time.
            
            - name: A user defined name for the reaction.
            - reactants: A dictionary with participant reactants as keys, and consumed per reaction as value.
            - products: A dictionary with reaction products as keys, and number formed per reaction as value.
            - rate: A parameter rate constant to be applied to the propensity of this reaction firing.
            - propensity_function: Can be used instead of rate in order to declare a custom propensity function in string format.
            """

            r1 = Reaction(
                    name="r1",
                    reactants={G3P: 1, inPD_R: 1}, 
                    products={G3P: 1, actPD:1, R:1},
                    rate=rate1
                )
            
            r2 = Reaction(
                    name="r2",
                    reactants={G3P: 1, inPK_R: 1}, 
                    products={G3P: 1, actPK:1, R:1},
                    rate=rate2
                )
            
            r3 = Reaction(
                    name="r3",
                    reactants={R:1, actPD: 1}, 
                    products={inPD_R: 1},
                    rate=rate3
                )
            
            r4 = Reaction(
                    name="r4",
                    reactants={R:1, actPK: 1}, 
                    products={inPK_R: 1},
                    rate=rate4
                )
            
            r5 = Reaction(
                    name="r5",
                    reactants={actPD: 1}, 
                    products={mD: 1, actPD:1},
                    rate=rate5
                )
            
            r6 = Reaction(
                    name="r6",
                    reactants={actPK: 1}, 
                    products={mK: 1, actPK:1},
                    rate=rate6
                )
            
            r7 = Reaction(
                    name="r7",
                    reactants={inPD_R: 1}, 
                    products={inPD_R: 1, mD: 1},
                    rate=rate7
                )
            
            r8 = Reaction(
                    name="r8",
                    reactants={inPK_R: 1}, 
                    products={inPK_R: 1, mK: 1},
                    rate=rate8
                )
                
            r9 = Reaction(
                    name="r9",
                    reactants={PR: 1}, 
                    products={PR: 1, mR: 1},
                    rate=rate9
                )
            
            r10 = Reaction(
                    name="r10",
                    reactants={mD: 1}, 
                    products={mD: 1, D: 1},
                    rate=rate10
                )
            
            r11 = Reaction(
                    name="r11",
                    reactants={mK: 1}, 
                    products={mK: 1, K: 1},
                    rate=rate11
                )
            r12 = Reaction(
                    name="r12",
                    reactants={mR: 1}, 
                    products={mR: 1, R: 1},
                    rate=rate11
                )
            r13 = Reaction(
                    name="r13",
                    reactants={}, 
                    products={G3P:1},
                    rate=rate13
                )
            r14 = Reaction(
                    name="r14",
                    reactants={}, 
                    products={Gly:1},
                    rate=rate14
                )        
            r15 = Reaction(
                    name="r15",
                    reactants={K: 1, Gly:1}, 
                    products={K: 1, G3P: 1},
                    rate=rate15
                )
            
            r16 = Reaction(
                    name="r16",
                    reactants={D: 1, G3P: 1}, 
                    products={D: 1},
                    rate=rate16
                )
            
            r17 = Reaction(
                    name="r17",
                    reactants={mD: 1}, 
                    products={},
                    rate=rate17
                )
            
            r18 = Reaction(
                    name="r18",
                    reactants={mK: 1}, 
                    products={},
                    rate=rate18
                )
            r19 = Reaction(
                    name="r19",
                    reactants={mR: 1}, 
                    products={},
                    rate=rate19
                )    
            
            r20 = Reaction(
                    name="r20",
                    reactants={D: 1}, 
                    products={},
                    rate=rate20
                )
            
            r21 = Reaction(
                    name="r21",
                    reactants={K: 1}, 
                    products={},
                    rate=rate21
                )
            
            r22 = Reaction(
                    name="r22",
                    reactants={R: 1}, 
                    products={},
                    rate=rate22
                )
            
            # Add the Reactions to the Model.
            self.add_reaction([r1, r2, r3, r4, r5, r6,
                               r7, r8, r9, r10, r11, r12,
                               r13, r14, r15, r16, r17, r18,
                               r19, r20, r21, r22])
            
            # Use NumPy to set the timespan of the Model.
            self.timespan(np.linspace(0, 1, 10))
            
            
# class AddValue:

#     def av(self):
#         model = MichaelisMenten(1,2,3,4,5,6,7,8,9)
#         results = model.run(solver=TauLeapingSolver)
#         return results
        

# for i in range(1):
#     model = MichaelisMenten(1,1,0,0,5,5,5,5,5,5,5,5,5,5)
#     results = model.run(solver=NumPySSASolver)
#     results.plot(title='testC.png',included_species_list=['mD','mK','mR'],save_png=True)
# plt.plot(results['time'],results['G3P'])
# print(np.mean(results['G3P']))
# b=AddValue()
# kk=b.av()

# kk.plot()





class Gillespie(Model):
     def __init__(self, pinpd, pinpk, pactpd, pactpk, ppr, pg3p, pmd, pmk, pmr, pd, pk, pr, pgly, pg3pr, ppar, parameter_values=None):

            # Intialize the Model with a name of your choosing.
            Model.__init__(self, name="Michaelis_Menten")
            
            """
            Parameters are constant values relevant to the system, such as reaction kinetic rates.
            
            - name: A user defined name for reference.
            - expression: Some constant value.
            """
            rate1 = Parameter(name="kbR", expression=0.1)
            rate2 = Parameter(name="kuR", expression=0.01)
            rate3 = Parameter(name="koffD", expression=0.1)
            rate4 = Parameter(name="konD", expression=0.1)
            rate5 = Parameter(name="koffK", expression=0.1)
            rate6 = Parameter(name="konK", expression=0.1)
            rate7 = Parameter(name="ktrsc", expression=0.1)
            rate8 = Parameter(name="ktrsc1", expression=0.1)
            rate9 = Parameter(name="ktrsc2", expression=0.01)
            rate10 = Parameter(name="ktrsc3", expression=0.01)
            rate11 = Parameter(name="konR", expression=0.5)
            rate12 = Parameter(name="ktrsl", expression=0.2)
            rate13 = Parameter(name="ktrsl1", expression=0.2)
            rate14 = Parameter(name="ktrsl2", expression=0.2)
            rate15 = Parameter(name="ktrsp_g3p", expression=0)
            rate16 = Parameter(name="ktrsp_gly", expression=0)
            rate17 = Parameter(name="kpf", expression=0.005)
            rate18 = Parameter(name="knf", expression=0.005)
            rate19 = Parameter(name="dmD", expression=0.069)
            rate20 = Parameter(name="dmK", expression=0.069)
            rate21 = Parameter(name="dmR", expression=0.069)
            rate22 = Parameter(name="dpD", expression=0.023)
            rate23 = Parameter(name="dpK", expression=0.023)
            rate24 = Parameter(name="dpR", expression=0.023)
            rate25 = Parameter(name="koffR", expression=0.1)
            rate26 = Parameter(name="ktr", expression=0.2)
            
            
            # Add the Parameters to the Model.
            self.add_parameter([rate1, rate2, rate3, rate4, rate5, rate6,
                                rate7, rate8, rate9, rate10, rate11, rate12,
                                rate13, rate14, rate15, rate16, rate17, rate18,
                                rate19, rate20, rate21, rate22, rate23, rate24, rate25, rate26])
            
            """
            Species can be anything that participates in or is produced by a reaction channel.
            
            - name: A user defined name for the species.
            - initial_value: A value/population count of species at start of simulation.
            """

            inPD_R = Species(name="inPD_R", initial_value=pinpd)
            inPK_R = Species(name="inPK_R", initial_value=pinpk)
            actPD = Species(name="actPD", initial_value=pactpd)
            actPK = Species(name="actPK", initial_value=pactpk)
            PR = Species(name="PR", initial_value=ppr)
            G3P = Species(name="G3P", initial_value=pg3p)
            mD = Species(name="mD", initial_value=pmd)
            mK = Species(name="mK", initial_value=pmk)
            mR = Species(name="mR", initial_value=pmr)
            D = Species(name="D", initial_value=pd)
            K = Species(name="K", initial_value=pk)
            R = Species(name="R", initial_value=pr)
            Gly = Species(name="Gly", initial_value=pgly)
            G3PR = Species(name="G3PR", initial_value=pg3pr)
            PAR = Species(name="PAR", initial_value=ppar)
            
            # Add the Species to the Model.
            self.add_species([inPD_R, inPK_R, actPD, actPK, PR, G3P, mD, mK, mR, D, K, R, Gly, G3PR, PAR])
            
            """
            Reactions are the reaction channels which cause the system to change over time.
            
            - name: A user defined name for the reaction.
            - reactants: A dictionary with participant reactants as keys, and consumed per reaction as value.
            - products: A dictionary with reaction products as keys, and number formed per reaction as value.
            - rate: A parameter rate constant to be applied to the propensity of this reaction firing.
            - propensity_function: Can be used instead of rate in order to declare a custom propensity function in string format.
            """

            r1 = Reaction(
                    name="r1",
                    reactants={G3P: 1, R: 1}, 
                    products={G3PR: 1},
                    rate=rate1
                )
            
            r2 = Reaction(
                    name="r2",
                    reactants={G3PR: 1}, 
                    products={G3P: 1, R:1},
                    rate=rate2
                )
            
            r3 = Reaction(
                    name="r3",
                    reactants={R:1, actPD: 1}, 
                    products={inPD_R: 1},
                    rate=rate3
                )
            r4 = Reaction(
                    name="r4",
                    reactants={inPD_R: 1}, 
                    products={R:1, actPD: 1},
                    rate=rate4
                )
            
            r5 = Reaction(
                    name="r5",
                    reactants={R:1, actPK: 1}, 
                    products={inPK_R: 1},
                    rate=rate5
                )
            r6 = Reaction(
                    name="r6",
                    reactants={inPK_R: 1}, 
                    products={R:1, actPK: 1},
                    rate=rate6
                )
            
            r7 = Reaction(
                    name="r7",
                    reactants={actPD: 1}, 
                    products={mD: 1, actPD:1},
                    rate=rate7
                )
            
            r8 = Reaction(
                    name="r8",
                    reactants={actPK: 1}, 
                    products={mK: 1, actPK:1},
                    rate=rate8
                )
            
            r9 = Reaction(
                    name="r9",
                    reactants={inPD_R: 1}, 
                    products={inPD_R: 1, mD: 1},
                    rate=rate9
                )
            
            r10 = Reaction(
                    name="r10",
                    reactants={inPK_R: 1}, 
                    products={inPK_R: 1, mK: 1},
                    rate=rate10
                )
                
            r11 = Reaction(
                    name="r11",
                    reactants={R:1, PR: 1}, 
                    products={PAR: 1},
                    rate=rate11
                )
            
            
            r12 = Reaction(
                    name="r12",
                    reactants={mD: 1}, 
                    products={mD: 1, D: 1},
                    rate=rate12
                )
            
            r13 = Reaction(
                    name="r13",
                    reactants={mK: 1}, 
                    products={mK: 1, K: 1},
                    rate=rate13
                )
            r14 = Reaction(
                    name="r14",
                    reactants={mR: 1}, 
                    products={mR: 1, R: 1},
                    rate=rate14
                )
            r15 = Reaction(
                    name="r15",
                    reactants={}, 
                    products={G3P:1},
                    rate=rate15
                )
            r16 = Reaction(
                    name="r16",
                    reactants={}, 
                    products={Gly:1},
                    rate=rate16
                )        
            r17 = Reaction(
                    name="r17",
                    reactants={K: 1, Gly:1}, 
                    products={K: 1, G3P: 1},
                    rate=rate17
                )
            
            r18 = Reaction(
                    name="r18",
                    reactants={D: 1, G3P: 1}, 
                    products={D: 1},
                    rate=rate18
                )
            
            r19 = Reaction(
                    name="r19",
                    reactants={mD: 1}, 
                    products={},
                    rate=rate19
                )
            
            r20 = Reaction(
                    name="r20",
                    reactants={mK: 1}, 
                    products={},
                    rate=rate20
                )
            r21 = Reaction(
                    name="r21",
                    reactants={mR: 1}, 
                    products={},
                    rate=rate21
                )    
            
            r22 = Reaction(
                    name="r22",
                    reactants={D: 1}, 
                    products={},
                    rate=rate22
                )
            
            r23 = Reaction(
                    name="r23",
                    reactants={K: 1}, 
                    products={},
                    rate=rate23
                )
            
            r24 = Reaction(
                    name="r24",
                    reactants={R: 1}, 
                    products={},
                    rate=rate24
                )
            r25 = Reaction(
                    name="r25",
                    reactants={PAR: 1}, 
                    products={R:1, PR: 1},
                    rate=rate25
                )
            
            r26 = Reaction(
                    name="r26",
                    reactants={PAR: 1}, 
                    products={PAR:1, mR: 1},
                    rate=rate26
                )
            
            # Add the Reactions to the Model.
            self.add_reaction([r1, r2, r3, r4, r5, r6,
                               r7, r8, r9, r10, r11, r12,
                               r13, r14, r15, r16, r17, r18,
                               r19, r20, r21, r22, r23, r24, r25, r26])
            
            # Use NumPy to set the timespan of the Model.
            self.timespan(np.linspace(0, 2, 20))
            




# import time
# start=time.time()
# # nd=np.array([])
# for i in range(1,5):
#     print(i)
#     model = Gillespie(1,1,0,0,0,5,5,5,5,5,5,5,5,5,1)
#     results = model.run(solver=NumPySSASolver)
#     # nd=np.append(nd,results['G3P'][-10:-1])
#     results.plot(title='testC.png',included_species_list=['mD'],save_png=True)


# print(time.time()-start)
# # model = Gillespie(1,1,0,0,0,5,5,5,5,5,5,5,5,5,1)
# results = model.run(solver=TauLeapingSolver)



# results.plot(title='plot1.png',included_species_list=['G3P'],save_png=True)

# results.plot(title='plot2.png',included_species_list=['mD','mK','mR'],save_png=True)

# results.plot(title='plot3.png',included_species_list=['D','K'],save_png=True)

# results.plot(included_species_list=['G3P'])

# plt.hist(results['G3P'], bins='auto')



# import seaborn as sns
# data = results['mD']
# sns.distplot(data, hist=False)











# plt.plot(results['time'],results['G3P'])
# plt.plot(results['time'],results['mD'])
# plt.savefig('plot.png')
# plt.close()
# print(np.mean(results['G3P']))
# b=AddValue()
# kk=b.av()

# kk.plot()
            
            
# Instantiate your Model.
# model = MichaelisMenten(1000,50,0,0)

"""
Run a stochastic simulation on the Model and store the results in the 'results' variable.
GillesPy2 will use the best solver for the Model if no solver is declared (see below).
"""

# results = model.run()


"""
These are the possible pure-Python solvers. These may be faster for small (tiny) problem sizes, 
but performance can drop off extremely quickly.
"""

# from gillespy2.solvers.numpy import (
# 	NumPySSASolver,
# 	ODESolver,
# 	TauLeapingSolver,
# 	TauHybridSolver
# )

"""
These are the C++ solvers. These can have slower startup time but will quickly outpace their
Python counterparts for larger problem sets.

Note: The usage of the C++ solver family requires dependencies 'GCC' and 'Make' to be installed
on your system. If neither of these can be found, GillesPy2 will select alternate solvers.
"""
# from gillespy2.solvers.cpp import (
# 	SSACSolver,
# 	ODECSolver,
# 	TauLeapingCSolver
# )

"""
The following code shows how to run the previously defined Model on a solver of your choosing.
Note: You may use the prior imports as a reference, but you do not need to import them all.

The following code shows two example usages -- one using the C++ SSA solver, the other 
using the Python Tau Leaping solver.
"""

# Import the C++ SSA Solver. This will compute highly accurate results, but will be slow for large 
# problem sets.
# from gillespy2.solvers.cpp import SSACSolver

# Run the previously defined and instantiate Model with an extra argument.
# results = model.run(solver=TauLeapingCSolver)

# If we instead want to use a Python solver, we can import and use it in the same way.
# from gillespy2.solvers.numpy import TauLeapingSolver

# results = model.run(solver=TauLeapingSolver)

"""
Plot the results of the simulation. 

There are a multitude of arguments that can be set to tweak the behavior and visuals of the plot. 
For now though, lets run it with default settings.
"""

# results.plot()
# plt.savefig('foo.pdf')
