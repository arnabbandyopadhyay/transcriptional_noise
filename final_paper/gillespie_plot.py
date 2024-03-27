# -*- coding: utf-8 -*-
"""
Created on Sun Nov 21 20:24:59 2021

@author: Arnab
"""

# Numpy is used to set the timespan of the Model.
import numpy as np
#import matplotlib.pyplot as plt
import random



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


class Gillespie_bistable(Model):
     def __init__(self, ppd, ppd_r, ppk, ppk_r,ppr, pg3p, pmd, pmr, pmk, pd, pr, pk, pg3pr, ppar, factor,mu,influx,
                  parameter_values=None):

            # Intialize the Model with a name of your choosing.
            Model.__init__(self, name="Michaelis_Menten")
            
            """
            Parameters are constant values relevant to the system, such as reaction kinetic rates.
            
            - name: A user defined name for reference.
            - expression: Some constant value.
            """
            scaling=1e9
            kk=float(factor)

            rate1 = Parameter(name="kbpd", expression=0.1)#2165.44314226246/scaling)
            rate2 = Parameter(name="kupd", expression=0.5)#10014.7094615103/scaling)
            rate3 = Parameter(name="kbpk", expression=0.1)#1588.06318455853/scaling)
            rate4 = Parameter(name="kupk", expression=0.5)#10477.3841801672/scaling)
            rate5 = Parameter(name="kbpr", expression=0)# 0.05 also gives interesting
            rate6 = Parameter(name="kupr", expression=0)
        
            rate7 = Parameter(name="ktrd", expression=1*influx)
            rate8 = Parameter(name="ktrk", expression=1*influx)
            rate9 = Parameter(name="ktrr", expression=1*influx/2)
            
            rate10 = Parameter(name="ktsd", expression=0.2)
            rate11 = Parameter(name="ktsk", expression=0.2)
            rate12 = Parameter(name="ktsr", expression=0.2)
            
            rate13 = Parameter(name="king3p", expression=0.5)#0.2, 0.01757360326999)
            rate14 = Parameter(name="kb", expression=0.1)#0.1#99.7805018542808/scaling) 
            rate15 = Parameter(name="ku", expression=0.001)
            
            rate16 = Parameter(name="kddg3p", expression=0.05)#*np.exp(-(mu-30)/10))#0.103423441900033)
                      
            rate17 = Parameter(name="kdmd", expression=(0.02+0.05))#*np.exp(-(mu-30)/50))
            rate18 = Parameter(name="kdmk", expression=(0.02+0.05))#*np.exp(-(mu-30)/50))
            rate19 = Parameter(name="kdmr", expression=(0.02+0.05))#*np.exp(-(mu-30)/50))            
            
            rate20 = Parameter(name="kdd", expression=(0.01+0.023))#*np.exp(-(mu-30)/50))
            rate21 = Parameter(name="kdk", expression=(0.01+0.023))#*np.exp(-(mu-30)/50))
            rate22 = Parameter(name="kdr", expression=(0.01+0.023))#*np.exp(-(mu-30)/50))     
            rate23 = Parameter(name="kdg3p", expression=(0.01+0.023))#*np.exp(-(mu-30)/50))
            
            rate24 = Parameter(name="kbasal", expression=(0.2*kk))#*np.exp(-(mu-30)/50))
           
            
            # Add the Parameters to the Model.
            self.add_parameter([rate1, rate2, rate3, rate4, rate5, rate6,
                                rate7, rate8, rate9, rate10, rate11, rate12,
                                rate13, rate14, rate15, rate16, rate17, rate18,
                                rate19, rate20, rate21, rate22, rate23, rate24])
            
            """
            Species can be anything that participates in or is produced by a reaction channel.
            
            - name: A user defined name for the species.
            - initial_value: A value/population count of species at start of simulation.
            """
            pD = Species(name="pD", initial_value=ppd)
            pD_R = Species(name="pD_R", initial_value=ppd_r)
            pK = Species(name="pK", initial_value=ppk)
            pK_R = Species(name="pK_R", initial_value=ppk_r)
            pR = Species(name="pR", initial_value=ppr)
            G3P = Species(name="G3P", initial_value=pg3p)
            mD = Species(name="mD", initial_value=pmd)
            mR = Species(name="mR", initial_value=pmr)
            mK = Species(name="mK", initial_value=pmk)
            D = Species(name="D", initial_value=pd)            
            R = Species(name="R", initial_value=pr)
            K = Species(name="K", initial_value=pk)
            G3PR = Species(name="G3PR", initial_value=pg3pr)
            paR = Species(name="paR", initial_value=ppar)
            
            
            
            # Add the Species to the Model.
            self.add_species([pD, pD_R, pK, pK_R, pR, G3P, mD, mR, mK, 
                              D, R, K, G3PR, paR])
            
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
                    reactants={R:1, pD: 1}, 
                    products={pD_R: 1},
                    rate=rate1
                )
            r2 = Reaction(
                    name="r2",
                    reactants={pD_R: 1}, 
                    products={R:1, pD: 1},
                    rate=rate2
                )
            
            r3 = Reaction(
                   name="r3",
                   reactants={R:1, pK: 1}, 
                   products={pK_R: 1},
                   rate=rate3
               )
            
            r4 = Reaction(
                    name="r4",
                    reactants={pK_R: 1}, 
                    products={R:1, pK: 1},
                    rate=rate4
                )
            r5 = Reaction(
                    name="r5",
                    reactants={paR: 1, R:1}, 
                    products={pR:1},
                    rate=rate5
                )
            r6 = Reaction(
                    name="r6",
                    reactants={pR:1}, 
                    products={paR: 1, R:1},
                    rate=rate6
                )
            r7 = Reaction(
                    name="r7",
                    reactants={pD: 1}, 
                    products={pD: 1, mD: 1},
                    rate=rate7
                )
            r8 = Reaction(
                    name="r8",
                    reactants={pK: 1}, 
                    products={pK: 1, mK: 1},
                    rate=rate8
                )
            r9 = Reaction(
                    name="r9",
                    reactants={pR: 1}, 
                    products={pR: 1, mR: 1},
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
                    rate=rate12
                )
            
            r13 = Reaction(
                    name="r13",
                    reactants={K:1}, 
                    products={G3P: 1, K:1},
                    rate=rate13
                )
            r14 = Reaction(
                    name="r14",
                    reactants={R: 1, G3P: 1}, 
                    products={G3PR:1},
                    rate=rate14
                )
            r15 = Reaction(
                    name="r15",
                    reactants={G3PR:1}, 
                    products={G3P: 1, R:1},
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
            r23 = Reaction(
                    name="r23",
                    reactants={G3P: 1}, 
                    products={},
                    rate=rate23
                )
            r24 = Reaction(
                    name="r24",
                    reactants={}, 
                    products={G3P:1},
                    rate=rate24
                )
            
            
    
            
         # Add the Reactions to the Model.
            self.add_reaction([r1, r2, r3, r4, r5, r6,
                               r7, r8, r9, r10, r11, r12,
                               r13, r14, r15, r16, r17, r18,
                               r19, r20, r21, r22, r23, r24])
            
            # Use NumPy to set the timespan of the Model.
            self.timespan(np.linspace(0, 200, 200))


import pandas as pd
cv = lambda x: np.std(x, ddof=1) / np.mean(x)
#import seaborn as sns
#import matplotlib.pyplot as plt
nd=np.empty((0,6),int)
par_var=[0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0]
#par_var=[ 0.7, 1.0]
for j in par_var:
   print(j)
   cv_d=np.empty((0,6),int)
   for i in range(1,500):
       
       if j<0.15:
       	model = Gillespie_bistable(1,0,1,0,1,5,5,5,5,5,5,5,5,0,0,30,j)
       # model = Gillespie_bistable(1,0,1,0,1,100,100,100,100,100,100,100,100,1,1)
       elif j<0.7:
       	model = Gillespie_bistable(1,0,1,0,1,5,5,5,5,5,5,5,5,0,0,30,j)
       else:
       	model = Gillespie_bistable(1,0,1,0,1,100,5,5,5,5,5,5,5,0,2,30,j)
     
       results = model.run(solver=NumPySSASolver)
       cv_d=np.append(cv_d,np.array([[j, results['mD'][-1],results['D'][-1],results['G3P'][-1],results['mR'][-1],results['pR'][-1]]]),axis=0)
       nd=np.append(nd,np.array([[j, results['mD'][-1],results['D'][-1],results['G3P'][-1],results['mR'][-1],results['pR'][-1]]]),axis=0)
   df=pd.DataFrame(cv_d)
   print(df.apply(cv))


df = pd.DataFrame(data=nd, columns=["param",'mD','D', 'G3P', 'mR','pR'])   
df.to_csv('test.csv',index=False)     

#df2=pd.melt(df, id_vars=['param'], value_vars=['mD', 'mR'])
#fig=sns.violinplot(data=df2, x='param', y='value',hue="variable", split=False, scale="count", cut=0,saturation=0.5)
#figure = fig.get_figure()    
#figure.savefig('violin.pdf', dpi=400)


