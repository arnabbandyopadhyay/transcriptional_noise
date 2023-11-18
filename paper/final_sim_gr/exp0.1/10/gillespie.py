# -*- coding: utf-8 -*-
"""
Created on Sun Nov 21 20:24:59 2021

@author: Arnab
"""

# Numpy is used to set the timespan of the Model.
import numpy as np
import matplotlib.pyplot as plt
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
            # rate1 = Parameter(name="kon_pd", expression=0.5)#2165.44314226246/scaling)
            # rate2 = Parameter(name="kon_pk", expression=0.05)#10014.7094615103/scaling)
            # rate3 = Parameter(name="koff_pr", expression=0.5)#1588.06318455853/scaling)
            # rate4 = Parameter(name="koff_pr2", expression=0.05)#10477.3841801672/scaling)
            
            # rate5 = Parameter(name="ktrsc1", expression=0.2*kk)
            # rate6 = Parameter(name="ktrsc2", expression=0.2*kk)
            # rate7 = Parameter(name="ktrsc3", expression=0.5)
            
            # rate8 = Parameter(name="ktrsl1", expression=1)
            # rate9 = Parameter(name="ktrsl2", expression=0.5)
            # rate10 = Parameter(name="ktrsl3", expression=0.5)
            
            # rate11 = Parameter(name="ktrsp_g3p", expression=0.15*kk)#0.01757360326999)
            # rate12 = Parameter(name="kon_g3pr", expression=0.005)#99.7805018542808/scaling)           
            # rate13 = Parameter(name="kon_g3pr2", expression=0.005*0)#0.103423441900033)
            # rate14 = Parameter(name="kon_g3pr3", expression=0.005*0)#99.7805018542808/scaling)           
            # rate15 = Parameter(name="kon_g3pd", expression=0.001)#0.103423441900033)
            # rate16 = Parameter(name="drna1", expression=0.069)
            # rate17 = Parameter(name="drna2", expression=0.069)
            # rate18 = Parameter(name="drna3", expression=0.069)            
            # rate19 = Parameter(name="dprotein1", expression=0.023)
            # rate20 = Parameter(name="dprotein2", expression=0.023)
            # rate21 = Parameter(name="dprotein3", expression=0.023)     
            # rate22 = Parameter(name="dg3p", expression=0.035)
            # rate23 = Parameter(name="kon_pr", expression=0.144)#2165.44314226246/scaling)
            # rate24 = Parameter(name="kof_pr", expression=0.5)#10014.7094615103/scaling)
            
            ############## current param
            
            # rate1 = Parameter(name="kbpd", expression=0.5)#2165.44314226246/scaling)
            # rate2 = Parameter(name="kupd", expression=0.05)#10014.7094615103/scaling)
            # rate3 = Parameter(name="kbpk", expression=0.5)#1588.06318455853/scaling)
            # rate4 = Parameter(name="kupk", expression=0.05)#10477.3841801672/scaling)
            # rate5 = Parameter(name="kbpr", expression=0.144)
            # rate6 = Parameter(name="kupr", expression=0.5)
        
            # rate7 = Parameter(name="ktrd", expression=0.2)
            # rate8 = Parameter(name="ktrk", expression=0.2)
            # rate9 = Parameter(name="ktrr", expression=0.05)
            
            # rate10 = Parameter(name="ktsd", expression=0.2)
            # rate11 = Parameter(name="ktsk", expression=0.2)
            # rate12 = Parameter(name="ktsr", expression=0.2)
            
            # rate13 = Parameter(name="king3p", expression=0.2)#0.2*kk*4, 0.01757360326999)
            # rate14 = Parameter(name="kb", expression=0.005)#99.7805018542808/scaling) 
            # rate15 = Parameter(name="ku", expression=0.001)
            
            # rate16 = Parameter(name="kddg3p", expression=0.5)#*np.exp(-(mu-30)/10))#0.103423441900033)
                      
            # rate17 = Parameter(name="kdmd", expression=5*0.069*kk)#np.exp(-(mu-30)/50))
            # rate18 = Parameter(name="kdmk", expression=5*0.069*kk)#np.exp(-(mu-30)/50))
            # rate19 = Parameter(name="kdmr", expression=5*0.069*kk)#np.exp(-(mu-30)/50))            
            
            # rate20 = Parameter(name="kdd", expression=0.023*kk)#np.exp(-(mu-30)/50))
            # rate21 = Parameter(name="kdk", expression=0.023*kk)#np.exp(-(mu-30)/50))
            # rate22 = Parameter(name="kdr", expression=0.023*kk)#np.exp(-(mu-30)/50))     
            # rate23 = Parameter(name="kdg3p", expression=0.023*kk)#np.exp(-(mu-30)/50))
    
            
            ######## param for xppaut bistability in R & D wrt ktrr glpD_simple.ode
            
            # rate1 = Parameter(name="kbpd", expression=0.05)#2165.44314226246/scaling)
            # rate2 = Parameter(name="kupd", expression=0.5)#10014.7094615103/scaling)
            # rate3 = Parameter(name="kbpk", expression=0.05)#1588.06318455853/scaling)
            # rate4 = Parameter(name="kupk", expression=0.5)#10477.3841801672/scaling)
            # rate5 = Parameter(name="kbpr", expression=0.05)
            # rate6 = Parameter(name="kupr", expression=0.5)
        
            # rate7 = Parameter(name="ktrd", expression=0.02*kk)
            # rate8 = Parameter(name="ktrk", expression=0.2*kk)
            # rate9 = Parameter(name="ktrr", expression=1.2)
            
            # rate10 = Parameter(name="ktsd", expression=0.2)
            # rate11 = Parameter(name="ktsk", expression=0.2)
            # rate12 = Parameter(name="ktsr", expression=0.2)
            
            # rate13 = Parameter(name="king3p", expression=0.5*kk)#0.01757360326999)
            # rate14 = Parameter(name="kb", expression=0.005)#99.7805018542808/scaling) 
            # rate15 = Parameter(name="ku", expression=0.001)
            
            # rate16 = Parameter(name="kddg3p", expression=0.05*np.exp(-(mu-30)/10))#0.103423441900033)
                      
            # rate17 = Parameter(name="kdmd", expression=0.069*np.exp(-(mu-30)/50))
            # rate18 = Parameter(name="kdmk", expression=0.069*np.exp(-(mu-30)/50))
            # rate19 = Parameter(name="kdmr", expression=0.069*np.exp(-(mu-30)/50))            
            
            # rate20 = Parameter(name="kdd", expression=0.023*np.exp(-(mu-30)/50))
            # rate21 = Parameter(name="kdk", expression=0.023*np.exp(-(mu-30)/50))
            # rate22 = Parameter(name="kdr", expression=0.023*np.exp(-(mu-30)/50))     
            # rate23 = Parameter(name="kdg3p", expression=0.023*np.exp(-(mu-30)/50))
            
            
            
            # ######## param for higher variability in D but lower in R (cell div 30)
            rate1 = Parameter(name="kbpd", expression=0.1)#2165.44314226246/scaling)
            rate2 = Parameter(name="kupd", expression=0.5)#10014.7094615103/scaling)
            rate3 = Parameter(name="kbpk", expression=0.1)#1588.06318455853/scaling)
            rate4 = Parameter(name="kupk", expression=0.5)#10477.3841801672/scaling)
            rate5 = Parameter(name="kbpr", expression=0)# 0.05 also gives interesting
            rate6 = Parameter(name="kupr", expression=0)
        
            rate7 = Parameter(name="ktrd", expression=1*influx)
            rate8 = Parameter(name="ktrk", expression=1*influx)
            rate9 = Parameter(name="ktrr", expression=0.05)
            
            rate10 = Parameter(name="ktsd", expression=0.2)
            rate11 = Parameter(name="ktsk", expression=0.2)
            rate12 = Parameter(name="ktsr", expression=0.2)
            
            rate13 = Parameter(name="king3p", expression=0.5)#0.2, 0.01757360326999)
            rate14 = Parameter(name="kb", expression=0.1)#0.1#99.7805018542808/scaling) 
            rate15 = Parameter(name="ku", expression=0.001)
            
            rate16 = Parameter(name="kddg3p", expression=0.05)#*np.exp(-(mu-30)/10))#0.103423441900033)
                      
            rate17 = Parameter(name="kdmd", expression=(0.02+0.05*kk))#*np.exp(-(mu-30)/50))
            rate18 = Parameter(name="kdmk", expression=(0.02+0.05*kk))#*np.exp(-(mu-30)/50))
            rate19 = Parameter(name="kdmr", expression=(0.02+0.05*kk))#*np.exp(-(mu-30)/50))            
            
            rate20 = Parameter(name="kdd", expression=(0.01+0.023*kk))#*np.exp(-(mu-30)/50))
            rate21 = Parameter(name="kdk", expression=(0.01+0.023*kk))#*np.exp(-(mu-30)/50))
            rate22 = Parameter(name="kdr", expression=(0.01+0.023*kk))#*np.exp(-(mu-30)/50))     
            rate23 = Parameter(name="kdg3p", expression=(0.01+0.023*kk))#*np.exp(-(mu-30)/50))
           
            
            # Add the Parameters to the Model.
            self.add_parameter([rate1, rate2, rate3, rate4, rate5, rate6,
                                rate7, rate8, rate9, rate10, rate11, rate12,
                                rate13, rate14, rate15, rate16, rate17, rate18,
                                rate19, rate20, rate21, rate22, rate23])
            
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
            
    
            
         # Add the Reactions to the Model.
            self.add_reaction([r1, r2, r3, r4, r5, r6,
                               r7, r8, r9, r10, r11, r12,
                               r13, r14, r15, r16, r17, r18,
                               r19, r20, r21, r22, r23])
            
            # Use NumPy to set the timespan of the Model.
            self.timespan(np.linspace(0, 200, 200))

class Gillespie_glycerol(Model):
     def __init__(self, ppd, ppd_r, ppk, ppk_r,ppr, pg3p, pmd, pmr, pmk, pd, pr, pk, pg3pr, ppar, factor,mu,
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
            # rate1 = Parameter(name="kon_pd", expression=0.5)#2165.44314226246/scaling)
            # rate2 = Parameter(name="kon_pk", expression=0.05)#10014.7094615103/scaling)
            # rate3 = Parameter(name="koff_pr", expression=0.5)#1588.06318455853/scaling)
            # rate4 = Parameter(name="koff_pr2", expression=0.05)#10477.3841801672/scaling)
            
            # rate5 = Parameter(name="ktrsc1", expression=0.2*kk)
            # rate6 = Parameter(name="ktrsc2", expression=0.2*kk)
            # rate7 = Parameter(name="ktrsc3", expression=0.5)
            
            # rate8 = Parameter(name="ktrsl1", expression=1)
            # rate9 = Parameter(name="ktrsl2", expression=0.5)
            # rate10 = Parameter(name="ktrsl3", expression=0.5)
            
            # rate11 = Parameter(name="ktrsp_g3p", expression=0.15*kk)#0.01757360326999)
            # rate12 = Parameter(name="kon_g3pr", expression=0.005)#99.7805018542808/scaling)           
            # rate13 = Parameter(name="kon_g3pr2", expression=0.005*0)#0.103423441900033)
            # rate14 = Parameter(name="kon_g3pr3", expression=0.005*0)#99.7805018542808/scaling)           
            # rate15 = Parameter(name="kon_g3pd", expression=0.001)#0.103423441900033)
            # rate16 = Parameter(name="drna1", expression=0.069)
            # rate17 = Parameter(name="drna2", expression=0.069)
            # rate18 = Parameter(name="drna3", expression=0.069)            
            # rate19 = Parameter(name="dprotein1", expression=0.023)
            # rate20 = Parameter(name="dprotein2", expression=0.023)
            # rate21 = Parameter(name="dprotein3", expression=0.023)     
            # rate22 = Parameter(name="dg3p", expression=0.035)
            # rate23 = Parameter(name="kon_pr", expression=0.144)#2165.44314226246/scaling)
            # rate24 = Parameter(name="kof_pr", expression=0.5)#10014.7094615103/scaling)
            
            ############## current param
            
            # rate1 = Parameter(name="kbpd", expression=0.5)#2165.44314226246/scaling)
            # rate2 = Parameter(name="kupd", expression=0.05)#10014.7094615103/scaling)
            # rate3 = Parameter(name="kbpk", expression=0.5)#1588.06318455853/scaling)
            # rate4 = Parameter(name="kupk", expression=0.05)#10477.3841801672/scaling)
            # rate5 = Parameter(name="kbpr", expression=0.144)
            # rate6 = Parameter(name="kupr", expression=0.5)
        
            # rate7 = Parameter(name="ktrd", expression=0.2)
            # rate8 = Parameter(name="ktrk", expression=0.2)
            # rate9 = Parameter(name="ktrr", expression=0.05)
            
            # rate10 = Parameter(name="ktsd", expression=0.2)
            # rate11 = Parameter(name="ktsk", expression=0.2)
            # rate12 = Parameter(name="ktsr", expression=0.2)
            
            # rate13 = Parameter(name="king3p", expression=0.2*kk*4)#0.2, 0.01757360326999)
            # rate14 = Parameter(name="kb", expression=0.005)#99.7805018542808/scaling) 
            # rate15 = Parameter(name="ku", expression=0.001)
            
            # rate16 = Parameter(name="kddg3p", expression=0.5)#*np.exp(-(mu-30)/10))#0.103423441900033)
                      
            # rate17 = Parameter(name="kdmd", expression=5*0.069)#np.exp(-(mu-30)/50))
            # rate18 = Parameter(name="kdmk", expression=5*0.069)#np.exp(-(mu-30)/50))
            # rate19 = Parameter(name="kdmr", expression=5*0.069)#np.exp(-(mu-30)/50))            
            
            # rate20 = Parameter(name="kdd", expression=0.023)#np.exp(-(mu-30)/50))
            # rate21 = Parameter(name="kdk", expression=0.023)#np.exp(-(mu-30)/50))
            # rate22 = Parameter(name="kdr", expression=0.023)#np.exp(-(mu-30)/50))     
            # rate23 = Parameter(name="kdg3p", expression=0.023)#np.exp(-(mu-30)/50))
    
            
            ######## param for xppaut bistability in R & D wrt ktrr glpD_simple.ode
            
            # rate1 = Parameter(name="kbpd", expression=0.05)#2165.44314226246/scaling)
            # rate2 = Parameter(name="kupd", expression=0.5)#10014.7094615103/scaling)
            # rate3 = Parameter(name="kbpk", expression=0.05)#1588.06318455853/scaling)
            # rate4 = Parameter(name="kupk", expression=0.5)#10477.3841801672/scaling)
            # rate5 = Parameter(name="kbpr", expression=0.05)
            # rate6 = Parameter(name="kupr", expression=0.5)
        
            # rate7 = Parameter(name="ktrd", expression=0.02*kk)
            # rate8 = Parameter(name="ktrk", expression=0.2*kk)
            # rate9 = Parameter(name="ktrr", expression=1.2)
            
            # rate10 = Parameter(name="ktsd", expression=0.2)
            # rate11 = Parameter(name="ktsk", expression=0.2)
            # rate12 = Parameter(name="ktsr", expression=0.2)
            
            # rate13 = Parameter(name="king3p", expression=0.5*kk)#0.01757360326999)
            # rate14 = Parameter(name="kb", expression=0.005)#99.7805018542808/scaling) 
            # rate15 = Parameter(name="ku", expression=0.001)
            
            # rate16 = Parameter(name="kddg3p", expression=0.05*np.exp(-(mu-30)/10))#0.103423441900033)
                      
            # rate17 = Parameter(name="kdmd", expression=0.069*np.exp(-(mu-30)/50))
            # rate18 = Parameter(name="kdmk", expression=0.069*np.exp(-(mu-30)/50))
            # rate19 = Parameter(name="kdmr", expression=0.069*np.exp(-(mu-30)/50))            
            
            # rate20 = Parameter(name="kdd", expression=0.023*np.exp(-(mu-30)/50))
            # rate21 = Parameter(name="kdk", expression=0.023*np.exp(-(mu-30)/50))
            # rate22 = Parameter(name="kdr", expression=0.023*np.exp(-(mu-30)/50))     
            # rate23 = Parameter(name="kdg3p", expression=0.023*np.exp(-(mu-30)/50))
            
            
            
            # ######## param for higher variability in D but lower in R (cell div 30)
            rate1 = Parameter(name="kbpd", expression=0.25)#2165.44314226246/scaling)
            rate2 = Parameter(name="kupd", expression=0.5)#10014.7094615103/scaling)
            rate3 = Parameter(name="kbpk", expression=0.25)#1588.06318455853/scaling)
            rate4 = Parameter(name="kupk", expression=0.5)#10477.3841801672/scaling)
            rate5 = Parameter(name="kbpr", expression=0.0)### 0.05 also gives interesting 
            rate6 = Parameter(name="kupr", expression=0.0)
        
            rate7 = Parameter(name="ktrd", expression=0.5)
            rate8 = Parameter(name="ktrk", expression=0.5)
            rate9 = Parameter(name="ktrr", expression=0.05)# 0.2-0.5
            
            rate10 = Parameter(name="ktsd", expression=0.2)
            rate11 = Parameter(name="ktsk", expression=0.2)
            rate12 = Parameter(name="ktsr", expression=0.2)
            
            rate13 = Parameter(name="king3p", expression=0.5*kk)#0.2, 0.01757360326999)
            rate14 = Parameter(name="kb", expression=0.1)#0.1#99.7805018542808/scaling) 
            rate15 = Parameter(name="ku", expression=0.001)
            
            rate16 = Parameter(name="kddg3p", expression=0.05)#*np.exp(-(mu-30)/10))#0.103423441900033)
                      
            rate17 = Parameter(name="kdmd", expression=0.069)#*np.exp(-(mu-30)/50))
            rate18 = Parameter(name="kdmk", expression=0.069)#*np.exp(-(mu-30)/50))
            rate19 = Parameter(name="kdmr", expression=0.069)#*np.exp(-(mu-30)/50))            
            
            rate20 = Parameter(name="kdd", expression=0.023)#*np.exp(-(mu-30)/50))
            rate21 = Parameter(name="kdk", expression=0.023)#*np.exp(-(mu-30)/50))
            rate22 = Parameter(name="kdr", expression=0.023)#*np.exp(-(mu-30)/50))     
            rate23 = Parameter(name="kdg3p", expression=0.023)#*np.exp(-(mu-30)/50))
           
            
            # Add the Parameters to the Model.
            self.add_parameter([rate1, rate2, rate3, rate4, rate5, rate6,
                                rate7, rate8, rate9, rate10, rate11, rate12,
                                rate13, rate14, rate15, rate16, rate17, rate18,
                                rate19, rate20, rate21, rate22, rate23])
            
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
            
    
            
         # Add the Reactions to the Model.
            self.add_reaction([r1, r2, r3, r4, r5, r6,
                               r7, r8, r9, r10, r11, r12,
                               r13, r14, r15, r16, r17, r18,
                               r19, r20, r21, r22, r23])
            
            # Use NumPy to set the timespan of the Model.
            self.timespan(np.linspace(0, 200, 200))

nd=np.empty((0,5),int)
nd2=np.empty((0,14),int)
for i in range(1,200):
    print(i)
    
    # model = Gillespie_bistable(1,0,1,0,1,100,100,100,100,100,100,100,100,1,1)
    model = Gillespie_bistable(1,0,1,0,1,5,5,5,5,5,5,5,5,0,1,30,0.4)
#    model = Gillespie_glycerol(1,0,1,0,1,1,1,1,1,1,1,1,1,0,0.5,30)
    results = model.run(solver=NumPySSASolver)
    nd=np.append(nd,np.array([[results['mD'][-1],results['D'][-1],results['G3P'][-1],results['R'][-1],results['pR'][-1]]]),axis=0)
    nd2=np.append(nd2,np.array([[results['pD'][-1],results['pD_R'][-1],results['pK'][-1],results['pK_R'][-1],
                                  results['pR'][-1], results['G3P'][-1], results['mD'][-1],
                                  results['mR'][-1], results['mK'][-1],results['D'][-1],
                                  results['R'][-1],results['K'][-1],results['G3PR'][-1],results['paR'][-1]]]),axis=0)
    # results.plot(included_species_list=['D'])
    
gly_dat=nd2
#import seaborn as sns
#data = nd[:,0] #results['mD']
#sns.distplot(data, hist=False, color='g')
## data = nd[:,3] #results['mD']
# sns.distplot(data, hist=False,color='r')
# model = Gillespie_bistable(1,0,1,0,1,1,1,1,1,1,1,1,1,1,1,60)
# results = model.run(solver=NumPySSASolver)
# results.plot(included_species_list=['D','R'])
# np.mean(data)
# np.std(data)
#plt.plot(data)

#import seaborn as sns
#nd=np.empty((0,6),int)
#par_var=[0.01, 0.03, 0.05, 0.08, 0.1, 0.15, 0.2]
#for j in par_var:
#    for i in range(1,1000):
#        print(i)
#        
#        # model = Gillespie_bistable(1,0,1,0,1,100,100,100,100,100,100,100,100,1,1)
#        model = Gillespie_bistable(1,0,1,0,1,5,5,5,5,5,5,5,5,0,1,30,j)
#      
#        results = model.run(solver=NumPySSASolver)
#        nd=np.append(nd,np.array([[j, results['mD'][-1],results['D'][-1],results['G3P'][-1],results['mR'][-1],results['pR'][-1]]]),axis=0)
#
#import pandas as pd
#df = pd.DataFrame(data=nd, columns=["param",'mD','D', 'G3P', 'mR','pR'])        
#
#df2=pd.melt(df, id_vars=['param'], value_vars=['mD', 'mR'])
#fig=sns.violinplot(data=df2, x='param', y='value',hue="variable", split=False, scale="count", cut=0)
#figure = fig.get_figure()    
#figure.savefig('violin.pdf', dpi=400)


# class Gillespie_bistable(Model):
#      def __init__(self, ppd, ppd_r, ppt, ppt_r,ppr, pg3p, pmd, pmr, pmt, pd, pr, pt, factor,
#                   parameter_values=None):

#             # Intialize the Model with a name of your choosing.
#             Model.__init__(self, name="Michaelis_Menten")
            
#             """
#             Parameters are constant values relevant to the system, such as reaction kinetic rates.
            
#             - name: A user defined name for reference.
#             - expression: Some constant value.
#             """
#             scaling=1e9
#             kk=float(factor)
#             rate1 = Parameter(name="kon_pdr", expression=0.144)#2165.44314226246/scaling)
#             rate2 = Parameter(name="koff_pdr", expression=0.0378)#10014.7094615103/scaling)
#             rate3 = Parameter(name="kon_ptr", expression=0.144)#1588.06318455853/scaling)
#             rate4 = Parameter(name="koff_ptr", expression=0.0378)#10477.3841801672/scaling)
#             rate5 = Parameter(name="ktrsc1", expression=0.5*kk)
#             rate6 = Parameter(name="ktrsc2", expression=0.5)
#             rate7 = Parameter(name="ktrsc3", expression=0.5/2)
#             rate8 = Parameter(name="ktrsl1", expression=0.3*kk)
#             rate9 = Parameter(name="ktrsl2", expression=0.3)
#             rate10 = Parameter(name="ktrsl3", expression=0.3)
#             rate11 = Parameter(name="ktrsp_g3p", expression=0.05*kk)#0.01757360326999)
#             rate12 = Parameter(name="kon_g3pr", expression=0.005)#99.7805018542808/scaling)           
#             rate13 = Parameter(name="kon_g3pd", expression=0.005*kk)#0.103423441900033)
#             rate14 = Parameter(name="drna1", expression=0.069)
#             rate15 = Parameter(name="drna2", expression=0.069)
#             rate16 = Parameter(name="drna3", expression=0.069)            
#             rate17 = Parameter(name="dprotein1", expression=0.023)
#             rate18 = Parameter(name="dprotein2", expression=0.023)
#             rate19 = Parameter(name="dprotein3", expression=0.023)     
#             rate20 = Parameter(name="dg3p", expression=0.035*kk)
            
            
            
#             # Add the Parameters to the Model.
#             self.add_parameter([rate1, rate2, rate3, rate4, rate5, rate6,
#                                 rate7, rate8, rate9, rate10, rate11, rate12,
#                                 rate13, rate14, rate15, rate16, rate17, rate18,
#                                 rate19, rate20])
            
#             """
#             Species can be anything that participates in or is produced by a reaction channel.
            
#             - name: A user defined name for the species.
#             - initial_value: A value/population count of species at start of simulation.
#             """
#             pD = Species(name="pD", initial_value=ppd)
#             pD_R = Species(name="pD_R", initial_value=ppd_r)
#             pT = Species(name="pT", initial_value=ppt)
#             pT_R = Species(name="pT_R", initial_value=ppt_r)
#             pR = Species(name="pR", initial_value=ppr)
#             G3P = Species(name="G3P", initial_value=pg3p)
#             mD = Species(name="mD", initial_value=pmd)
#             mR = Species(name="mR", initial_value=pmr)
#             mT = Species(name="mT", initial_value=pmt)
#             D = Species(name="D", initial_value=pd)            
#             R = Species(name="R", initial_value=pr)
#             T = Species(name="T", initial_value=pt)
            
            
            
#             # Add the Species to the Model.
#             self.add_species([pD, pD_R, pT, pT_R, pR, G3P, mD, mR, mT, 
#                               D, R, T])
            
#             """
#             Reactions are the reaction channels which cause the system to change over time.
            
#             - name: A user defined name for the reaction.
#             - reactants: A dictionary with participant reactants as keys, and consumed per reaction as value.
#             - products: A dictionary with reaction products as keys, and number formed per reaction as value.
#             - rate: A parameter rate constant to be applied to the propensity of this reaction firing.
#             - propensity_function: Can be used instead of rate in order to declare a custom propensity function in string format.
            
#             """
            
#             r1 = Reaction(
#                     name="r1",
#                     reactants={R:1, pD: 1}, 
#                     products={pD_R: 1},
#                     rate=rate1
#                 )
#             r2 = Reaction(
#                     name="r2",
#                     reactants={pD_R: 1}, 
#                     products={R:1, pD: 1},
#                     rate=rate2
#                 )
            
            
#             r3 = Reaction(
#                     name="r3",
#                     reactants={R:1, pT: 1}, 
#                     products={pT_R: 1},
#                     rate=rate3
#                 )
            
#             r4 = Reaction(
#                     name="r4",
#                     reactants={pT_R: 1}, 
#                     products={R:1, pT: 1},
#                     rate=rate4
#                 )
            
#             r5 = Reaction(
#                     name="r5",
#                     reactants={pD: 1}, 
#                     products={pD: 1, mD: 1},
#                     rate=rate5
#                 )
            
            
#             r6 = Reaction(
#                     name="r6",
#                     reactants={pT: 1}, 
#                     products={pT: 1, mT: 1},
#                     rate=rate6
#                 )
            
            
#             r7 = Reaction(
#                     name="r7",
#                     reactants={pR: 1}, 
#                     products={pR: 1, mR: 1},
#                     rate=rate7
#                 )
            
#             r8 = Reaction(
#                     name="r8",
#                     reactants={mD: 1}, 
#                     products={mD: 1, D: 1},
#                     rate=rate8
#                 )
            
            
            
#             r9 = Reaction(
#                     name="r9",
#                     reactants={mT: 1}, 
#                     products={mT: 1, T: 1},
#                     rate=rate9
#                 )
            
#             r10 = Reaction(
#                     name="r10",
#                     reactants={mR: 1}, 
#                     products={mR: 1, R: 1},
#                     rate=rate10
#                 )
            
#             r11 = Reaction(
#                     name="r11",
#                     reactants={T:1}, 
#                     products={G3P: 1, T:1},
#                     rate=rate11
#                 )
            
#             r12 = Reaction(
#                     name="r12",
#                     reactants={R: 1, G3P: 1}, 
#                     products={},
#                     rate=rate12
#                 )

      
#             r13 = Reaction(
#                     name="r13",
#                     reactants={D: 1, G3P: 1}, 
#                     products={D: 1},
#                     rate=rate13
#                 )
            
#             r14 = Reaction(
#                     name="r14",
#                     reactants={mD: 1}, 
#                     products={},
#                     rate=rate14
#                 )
#             r15 = Reaction(
#                     name="r15",
#                     reactants={mT: 1}, 
#                     products={},
#                     rate=rate15
#                 )  
 
#             r16 = Reaction(
#                     name="r16",
#                     reactants={mR: 1}, 
#                     products={},
#                     rate=rate16
#                 )                  
            
#             r17 = Reaction(
#                     name="r17",
#                     reactants={D: 1}, 
#                     products={},
#                     rate=rate17
#                 )
#             r18 = Reaction(
#                     name="r18",
#                     reactants={T: 1}, 
#                     products={},
#                     rate=rate18
#                 )

#             r19 = Reaction(
#                     name="r19",
#                     reactants={R: 1}, 
#                     products={},
#                     rate=rate19
#                 )

            
#             r20 = Reaction(
#                     name="r20",
#                     reactants={G3P: 1}, 
#                     products={},
#                     rate=rate20
#                 )
            
#          # Add the Reactions to the Model.
#             self.add_reaction([r1, r2, r3, r4, r5, r6,
#                                r7, r8, r9, r10, r11, r12,
#                                r13, r14, r15, r16, r17, r18,
#                                r19, r20])
            
#             # Use NumPy to set the timespan of the Model.
#             self.timespan(np.linspace(0, 10, 100))


# # import time
# # start=time.time()
# nd=np.empty((0,3),int)
# for i in range(1,300):
#     print(i)
    
#     # model = Gillespie_bistable(1,0,1,0,1,100,100,100,100,100,100,100)
#     model = Gillespie_bistable(1,0,1,0,1,9,9,10,9,9,30,9,20)
#     results = model.run(solver=NumPySSASolver)
#     nd=np.append(nd,np.array([[results['mD'][-1],results['D'][-1],results['G3P'][-1]]]),axis=0)
#     # results.plot(included_species_list=['D'])
# import seaborn as sns
# data = nd[:,0] #results['mD']
# sns.distplot(data, hist=False)

# results.plot(included_species_list=['mD'])

# print(time.time()-start)
# # model = Gillespie(1,1,0,0,0,5,5,5,5,5,5,5,5,5,1)
# results = model.run(solver=TauLeapingSolver)



# results.plot(title='plot1.png',included_species_list=['G3P'],save_png=True)

# results.plot(title='plot2.png',included_species_list=['mD','mK','mR'],save_png=True)

# results.plot(title='plot3.png',included_species_list=['D','K'],save_png=True)

# results.plot(included_species_list=['G3P'])

# plt.hist(results['G3P'], bins='auto')




# sns.distplot(results['mD'], hist=False)


# class Gillespie_bistable(Model):
#      def __init__(self, ppd, ppd_r, pfk, pfk_r, ppt, ppt_r,
#                   ppr, pg3p, pmd, pmk, pmr, pmf, pmt, pd, pk, pr, pf, pt, 
#                   pgly, pglyk, pg3pr, pg3pd, ita, parameter_values=None):

#             # Intialize the Model with a name of your choosing.
#             Model.__init__(self, name="Michaelis_Menten")
            
#             """
#             Parameters are constant values relevant to the system, such as reaction kinetic rates.
            
#             - name: A user defined name for reference.
#             - expression: Some constant value.
#             """
#             rate1 = Parameter(name="kon_pdr", expression=12.0999)
#             rate2 = Parameter(name="koff_pdr", expression=86.7934)
#             rate3 = Parameter(name="kon_pfkr", expression=11.0653)
#             rate4 = Parameter(name="koff_pfkr", expression=79.1819)
#             rate5 = Parameter(name="kon_ptr", expression=7.10524)
#             rate6 = Parameter(name="koff_ptr", expression=37.8405)
#             rate7 = Parameter(name="ktrsc", expression=0.5)
#             rate8 = Parameter(name="ktrsc1", expression=0.5)
#             rate9 = Parameter(name="ktrsc2", expression=0.5)
#             rate10 = Parameter(name="ktrsc3", expression=0.5/20)
#             rate11 = Parameter(name="ktrsc4", expression=0.5/20)
#             rate12 = Parameter(name="ktrsc5", expression=0.5/20)
#             rate13 = Parameter(name="ktrsc6", expression=0.5/2)           
#             rate14 = Parameter(name="ktrsl", expression=0.3)
#             rate15 = Parameter(name="ktrsl1", expression=0.3)
#             rate16 = Parameter(name="ktrsl2", expression=0.3)
#             rate17 = Parameter(name="ktrsl3", expression=0.3)
#             rate18 = Parameter(name="ktrsl4", expression=0.3)           
#             rate19 = Parameter(name="ktrsp_g3p", expression=0.0100652)
#             rate20 = Parameter(name="ktrsp_gly", expression=0.015)
#             rate21 = Parameter(name="kon_g3pr", expression=11.4709)
#             rate22 = Parameter(name="koff_g3pr", expression=0.360667)
#             rate23 = Parameter(name="kon_glyk", expression=37.1818)
#             rate24 = Parameter(name="koff_glyk", expression=0.794284)
#             rate25 = Parameter(name="kc_glyk", expression=7.5)
#             rate26 = Parameter(name="kon_g3pd", expression=0.01)
#             rate27 = Parameter(name="koff_g3pd", expression=0.278724)
#             rate28 = Parameter(name="kc_g3pd", expression=82.4)
#             rate29 = Parameter(name="drna", expression=0.069)
#             rate30 = Parameter(name="drna1", expression=0.069)
#             rate31 = Parameter(name="drna2", expression=0.069)
#             rate32 = Parameter(name="drna3", expression=0.069)
#             rate33 = Parameter(name="drna4", expression=0.069)
#             rate34 = Parameter(name="dprotein", expression=0.023)
#             rate35 = Parameter(name="dprotein1", expression=0.023)
#             rate36 = Parameter(name="dprotein2", expression=0.023)
#             rate37 = Parameter(name="dprotein3", expression=0.023)
#             rate38 = Parameter(name="dprotein4", expression=0.023)
#             rate39 = Parameter(name="dg3p", expression=0.035)
#             rate40 = Parameter(name="dg3pr", expression=0.035)
            
            
#             # Add the Parameters to the Model.
#             self.add_parameter([rate1, rate2, rate3, rate4, rate5, rate6,
#                                 rate7, rate8, rate9, rate10, rate11, rate12,
#                                 rate13, rate14, rate15, rate16, rate17, rate18,
#                                 rate19, rate20, rate21, rate22, rate23, rate24, rate25, rate26,
#                                 rate27, rate28, rate29, rate30, rate31, rate32,
#                                 rate33, rate34, rate35, rate36, rate37, rate38, rate39, rate40])
            
#             """
#             Species can be anything that participates in or is produced by a reaction channel.
            
#             - name: A user defined name for the species.
#             - initial_value: A value/population count of species at start of simulation.
#             """
#             pD = Species(name="pD", initial_value=ppd)
#             pD_R = Species(name="pD_R", initial_value=ppd_r)
#             pFK = Species(name="pFK", initial_value=pfk)
#             pFK_R = Species(name="pFK_R", initial_value=pfk_r)
#             pT = Species(name="pT", initial_value=ppt)
#             pT_R = Species(name="pT_R", initial_value=ppt_r)
#             pR = Species(name="pR", initial_value=ppr)
#             G3P = Species(name="G3P", initial_value=pg3p)
#             mD = Species(name="mD", initial_value=pmd)
#             mK = Species(name="mK", initial_value=pmk)
#             mR = Species(name="mR", initial_value=pmr)
#             mF = Species(name="mF", initial_value=pmf)
#             mT = Species(name="mT", initial_value=pmt)
#             D = Species(name="D", initial_value=pd)
#             K = Species(name="K", initial_value=pk)
#             R = Species(name="R", initial_value=pr)
#             F = Species(name="F", initial_value=pf)
#             T = Species(name="T", initial_value=pt)
#             Gly = Species(name="Gly", initial_value=pgly)
#             GlyK = Species(name="GlyK", initial_value=pglyk)
#             G3PR = Species(name="G3PR", initial_value=pg3pr)
#             G3PD = Species(name="G3PD", initial_value=pg3pd)
            
            
#             # Add the Species to the Model.
#             self.add_species([pD, pD_R, pFK, pFK_R, pT, pT_R, pR, G3P, mD, mK, mR, mF, mT, 
#                               D, K, R, F, T, Gly, GlyK, G3PR, G3PD])
            
#             """
#             Reactions are the reaction channels which cause the system to change over time.
            
#             - name: A user defined name for the reaction.
#             - reactants: A dictionary with participant reactants as keys, and consumed per reaction as value.
#             - products: A dictionary with reaction products as keys, and number formed per reaction as value.
#             - rate: A parameter rate constant to be applied to the propensity of this reaction firing.
#             - propensity_function: Can be used instead of rate in order to declare a custom propensity function in string format.
            
#             """
            
#             r1 = Reaction(
#                     name="r1",
#                     reactants={R:1, pD: 1}, 
#                     products={pD_R: 1},
#                     rate=rate1
#                 )
#             r2 = Reaction(
#                     name="r2",
#                     reactants={pD_R: 1}, 
#                     products={R:1, pD: 1},
#                     rate=rate2
#                 )
            
#             r3 = Reaction(
#                     name="r3",
#                     reactants={R:1, pFK: 1}, 
#                     products={pFK_R: 1},
#                     rate=rate3
#                 )
            
#             r4 = Reaction(
#                     name="r4",
#                     reactants={pFK_R: 1}, 
#                     products={R:1, pFK: 1},
#                     rate=rate4
#                 )
            
#             r5 = Reaction(
#                     name="r5",
#                     reactants={R:1, pT: 1}, 
#                     products={pT_R: 1},
#                     rate=rate5
#                 )
            
#             r6 = Reaction(
#                     name="r6",
#                     reactants={pT_R: 1}, 
#                     products={R:1, pT: 1},
#                     rate=rate6
#                 )
            
#             r7 = Reaction(
#                     name="r7",
#                     reactants={pD: 1}, 
#                     products={pD: 1, mD: 1},
#                     rate=rate7
#                 )
            
#             r8 = Reaction(
#                     name="r8",
#                     reactants={pFK: 1}, 
#                     products={pFK: 1, mF: 1, mK: 1},
#                     rate=rate8
#                 )
            
#             r9 = Reaction(
#                     name="r9",
#                     reactants={pT: 1}, 
#                     products={pT: 1, mT: 1},
#                     rate=rate9
#                 )
            
#             r10 = Reaction(
#                     name="r10",
#                     reactants={pD_R: 1}, 
#                     products={pD_R: 1, mD: 1},
#                     rate=rate10
#                 )
            
#             r11 = Reaction(
#                     name="r11",
#                     reactants={pFK_R: 1}, 
#                     products={pFK_R: 1, mF: 1, mK: 1},
#                     rate=rate11
#                 )
            
#             r12 = Reaction(
#                     name="r12",
#                     reactants={pT_R: 1}, 
#                     products={pT_R: 1, mT: 1},
#                     rate=rate12
#                 )
            
#             r13 = Reaction(
#                     name="r13",
#                     reactants={pR: 1}, 
#                     products={pR: 1, mR: 1},
#                     rate=rate13
#                 )
            
#             r14 = Reaction(
#                     name="r14",
#                     reactants={mD: 1}, 
#                     products={mD: 1, D: 1},
#                     rate=rate14
#                 )
            
#             r15 = Reaction(
#                     name="r15",
#                     reactants={mF: 1}, 
#                     products={mF: 1, F: 1},
#                     rate=rate15
#                 )
            
#             r16 = Reaction(
#                     name="r16",
#                     reactants={mK: 1}, 
#                     products={mK: 1, K: 1},
#                     rate=rate16
#                 )
            
#             r17 = Reaction(
#                     name="r17",
#                     reactants={mT: 1}, 
#                     products={mT: 1, T: 1},
#                     rate=rate17
#                 )
            
#             r18 = Reaction(
#                     name="r18",
#                     reactants={mR: 1}, 
#                     products={mR: 1, R: 1},
#                     rate=rate18
#                 )
            
#             r19 = Reaction(
#                     name="r19",
#                     reactants={T:1}, 
#                     products={T:1, G3P:1},
#                     rate=rate19
#                 )
#             r20 = Reaction(
#                     name="r20",
#                     reactants={F:1}, 
#                     products={F:1, Gly:1},
#                     rate=rate20
#                 )   
            
#             r21 = Reaction(
#                     name="r21",
#                     reactants={G3P: 1, R: 1}, 
#                     products={G3PR: 1},
#                     rate=rate21
#                 )
            
#             r22 = Reaction(
#                     name="r22",
#                     reactants={G3PR: 1}, 
#                     products={G3P: 1, R:1},
#                     rate=rate22
#                 )
            
#             r23 = Reaction(
#                     name="r23",
#                     reactants={Gly: 1, K: 1}, 
#                     products={GlyK: 1},
#                     rate=rate23
#                 )
            
#             r24 = Reaction(
#                     name="r24",
#                     reactants={GlyK: 1}, 
#                     products={Gly: 1, K: 1},
#                     rate=rate24
#                 )
            
#             r25 = Reaction(
#                     name="r25",
#                     reactants={GlyK: 1}, 
#                     products={G3P: 1, K: 1},
#                     rate=rate25
#                 )
            
#             r26 = Reaction(
#                     name="r26",
#                     reactants={D: 1, G3P: 1}, 
#                     products={G3PD: 1},
#                     rate=rate26
#                 )
#             r27 = Reaction(
#                     name="r27",
#                     reactants={G3PD: 1}, 
#                     products={D: 1, G3P: 1},
#                     rate=rate27
#                 )
#             r28 = Reaction(
#                     name="r28",
#                     reactants={G3PD: 1}, 
#                     products={D: 1},
#                     rate=rate28
#                 )
            
#             r29 = Reaction(
#                     name="r29",
#                     reactants={mD: 1}, 
#                     products={},
#                     rate=rate29
#                 )
            
#             r30 = Reaction(
#                     name="r30",
#                     reactants={mK: 1}, 
#                     products={},
#                     rate=rate30
#                 )
#             r31 = Reaction(
#                     name="r31",
#                     reactants={mR: 1}, 
#                     products={},
#                     rate=rate31
#                 )    
            
#             r32 = Reaction(
#                     name="r32",
#                     reactants={mF: 1}, 
#                     products={},
#                     rate=rate32
#                 )  
#             r33 = Reaction(
#                     name="r33",
#                     reactants={mT: 1}, 
#                     products={},
#                     rate=rate33
#                 )    
            
#             r34 = Reaction(
#                     name="r34",
#                     reactants={D: 1}, 
#                     products={},
#                     rate=rate34
#                 )
            
#             r35 = Reaction(
#                     name="r35",
#                     reactants={K: 1}, 
#                     products={},
#                     rate=rate35
#                 )
            
#             r36 = Reaction(
#                     name="r36",
#                     reactants={R: 1}, 
#                     products={},
#                     rate=rate36
#                 )
#             r37 = Reaction(
#                     name="r37",
#                     reactants={F: 1}, 
#                     products={},
#                     rate=rate37
#                 )
#             r38 = Reaction(
#                     name="r38",
#                     reactants={T: 1}, 
#                     products={},
#                     rate=rate38
#                 )
#             r39 = Reaction(
#                     name="r39",
#                     reactants={G3P: 1}, 
#                     products={},
#                     rate=rate39
#                 )
#             r40 = Reaction(
#                     name="r40",
#                     reactants={G3PR: 1}, 
#                     products={},
#                     rate=rate40
#                 )
#          # Add the Reactions to the Model.
#             self.add_reaction([r1, r2, r3, r4, r5, r6,
#                                r7, r8, r9, r10, r11, r12,
#                                r13, r14, r15, r16, r17, r18,
#                                r19, r20, r21, r22, r23, r24, r25, r26,
#                                r27,r28,r29,r30,r31,r32,r33,r34,r35,
#                                r36,r37,r38,r39,r40])
            
#             # Use NumPy to set the timespan of the Model.
#             self.timespan(np.linspace(0, 5, 50))











# import time
# start=time.time()
# nd=np.empty((0,3),int)
# for i in range(1,100):
#     print(i)
#     # model = Gillespie_bistable(1,0,1,0,1,0,1,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5)
#     model = Gillespie_bistable(1,0,1,0,1,0,1,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10)
#     results = model.run(solver=NumPySSASolver)
#     nd=np.append(nd,np.array([[results['mD'][-1],results['D'][-1],results['G3P'][-1]]]),axis=0)
#     # results.plot(title='testC.png',included_species_list=['G3P'],save_png=True)

# results.plot(included_species_list=['mD'])

# print(time.time()-start)
# # model = Gillespie(1,1,0,0,0,5,5,5,5,5,5,5,5,5,1)
# results = model.run(solver=TauLeapingSolver)



# results.plot(title='plot1.png',included_species_list=['G3P'],save_png=True)

# results.plot(title='plot2.png',included_species_list=['mD','mK','mR'],save_png=True)

# results.plot(title='plot3.png',included_species_list=['D','K'],save_png=True)

# results.plot(included_species_list=['G3P'])

# plt.hist(results['G3P'], bins='auto')



# import seaborn as sns
# data = nd[:,0] #results['mD']
# sns.distplot(data, hist=False)

# sns.distplot(results['mD'], hist=False)





















class Gillespie_old(Model):
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
            rate16 = Parameter(name="ktrsp_gly", expression=1)
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
            self.timespan(np.linspace(0, 1000, 1000))
            

# import time
# start=time.time()
# nd=np.empty((0,3),int)
# for i in range(1,500):
#     print(i)
#     model = Gillespie_old(1,1,0,0,0,5,5,5,5,5,5,5,5,5,1)
#     results = model.run(solver=NumPySSASolver)
#     nd=np.append(nd,np.array([[results['mD'][-1],results['D'][-1],results['G3P'][-1]]]),axis=0)
#     # results.plot(title='testC.png',included_species_list=['G3P'],save_png=True)

# # import seaborn as sns
# data = nd[:,1] #results['mD']
# sns.distplot(data, hist=False)



class Gillespie(Model):
     def __init__(self, pinpd, pinpk, pactpd, pactpk, ppr, pg3p, pmd, 
                  pmk, pmr, pd, pk, pr, pgly, pg3pr, ppar, ita, parameter_values=None):

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
            rate7 = Parameter(name="ktrsc", expression=1)
            rate8 = Parameter(name="ktrsc1", expression=1)
            rate9 = Parameter(name="ktrsc2", expression=0.1+0.1*ita)
            rate10 = Parameter(name="ktrsc3", expression=0.1)
            rate11 = Parameter(name="konR", expression=0.5)
            rate12 = Parameter(name="ktrsl", expression=0.2)
            rate13 = Parameter(name="ktrsl1", expression=0.2)
            rate14 = Parameter(name="ktrsl2", expression=0.2)
            rate15 = Parameter(name="ktrsp_g3p", expression=1+ita)
            rate16 = Parameter(name="ktrsp_gly", expression=1)
            rate17 = Parameter(name="kpf", expression=0.005)
            rate18 = Parameter(name="knf", expression=0.0001)
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
            self.timespan(np.linspace(0, 1000, 1000))
            




# import time
# start=time.time()
# nd=np.empty((0,3),int)
# for i in range(1,500):
#     print(i)
#     model = Gillespie(1,1,0,0,0,5,5,5,5,5,5,5,5,5,1,1)
#     # model = Gillespie(1,1,0,0,0,50,50,50,50,50,50,50,50,50,1,1)
#     results = model.run(solver=NumPySSASolver)
#     nd=np.append(nd,np.array([[results['mD'][-1],results['D'][-1],results['G3P'][-1]]]),axis=0)
#     # results.plot(title='testC.png',included_species_list=['G3P'],save_png=True)


# print(time.time()-start)
# # model = Gillespie(1,1,0,0,0,5,5,5,5,5,5,5,5,5,1)
# results = model.run(solver=TauLeapingSolver)



# results.plot(title='plot1.png',included_species_list=['G3P'],save_png=True)

# results.plot(title='plot2.png',included_species_list=['mD','mK','mR'],save_png=True)

# results.plot(title='plot3.png',included_species_list=['D','K'],save_png=True)

# results.plot(included_species_list=['mD'])

# plt.hist(results['G3P'], bins='auto')



# import seaborn as sns
# data = nd[:,2] #results['mD']
# sns.distplot(data, hist=False)

# sns.distplot(results['mD'], hist=False)











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
