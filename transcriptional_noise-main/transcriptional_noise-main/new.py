# -*- coding: utf-8 -*-
"""
Created on Sun Nov 21 10:58:43 2021

@author: Arnab
"""

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib import animation
import numpy as np
import random
from celluloid import Camera
import numpy as np
import gillespie as gp
from gillespy2.solvers.cpp import SSACSolver

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

mu, sigma = 120, 3 

class Cell:
    
    def __init__(self, birthtime, **kwargs):
        
        
        self.birthtime=birthtime
        
         
        self.nextdiv=self.birthtime+np.random.normal(mu, sigma)#np.random.normal(mu, sigma)
        
        if 'pinpd' in kwargs: self.pinpd = kwargs['pinpd']
        else: self.pinpd = 1
        
        if 'pinpk' in kwargs: self.pinpk = kwargs['pinpk']
        else: self.pinpk = 1
        
        if 'pactpd' in kwargs: self.pactpd = kwargs['pactpd']
        else: self.pactpd = 0
        
        if 'pactpk' in kwargs: self.pactpk = kwargs['pactpk']
        else: self.pactpk = 0
        
        
        if 'pmd' in kwargs: self.pmd = kwargs['pmd']
        else: self.pmd = np.random.randint(5,15)
        
        if 'pmk' in kwargs: self.pmk = kwargs['pmk']
        else: self.pmk = np.random.randint(5,15)
        
        if 'pd' in kwargs: self.pd = kwargs['pd']
        else: self.pd = np.random.randint(5,15)
        
        if 'pk' in kwargs: self.pk = kwargs['pk']
        else: self.pk = np.random.randint(5,15)
        
        if 'pg3p' in kwargs: self.pg3p = kwargs['pg3p']
        else: self.pg3p = np.random.randint(5,15)
        
        # self.generation = 0
        
        if 'xpos' in kwargs: self.xpos = kwargs['xpos']
        else: self.xpos = round(np.random.random(),2)*box_width
        
        if 'ypos' in kwargs: self.ypos = kwargs['ypos']
        else: self.ypos =round(np.random.random(),2)*box_width
        
        if 'xvel' in kwargs: self.xvel = kwargs['xvel']
        else: self.xvel = 2*(np.random.random()-0.5)/10
        
        if 'yvel' in kwargs: self.yvel = kwargs['yvel']
        else: self.yvel =2*(np.random.random()-0.5)/10
        
        
        
            
        
    @classmethod  
    def clone(cls, b,md,mk,d,k,g3p):
        return cls(birthtime=tnow, nextdiv=b.tnow+np.random.normal(mu, sigma),
                   pinpd = 1, pinpk = 1, pactpd = 0, pactpk = 0,
                   pmd=md,pmk=mk,pd=d,pk=k,pg3p=g3p,
                   xpos = b.xpos,
                   ypos = b.ypos,
                   xvel=2*(np.random.random()-0.5)/2,
                   yvel=2*(np.random.random()-0.5)/2)
    def gen(self):
        self.generation += 1

    def proliferate(self):


        nl=[Cell.clone(self),Cell.clone(self)]

        return nl
    
    def ssa(self):
        model = gp.MichaelisMenten(self.pinpd,self.pinpk,self.pactpd,
                                   self.pactpk,self.pmd, self.pmk, 
                                   self.pd, self.pk, self.pg3p)
        results = model.run(TauLeapingSolver)
        self.pinpd=results['inPD'][-1]
        self.pinpk=results['inPK'][-1]
        self.pactpd=results['actPD'][-1]
        self.pactpk=results['actPK'][-1]
        self.pmd=results['mD'][-1]
        self.pmk=results['mK'][-1]
        self.pd=results['D'][-1]
        self.pk=results['K'][-1]
        self.pg3p=results['G3P'][-1]
        


fig = plt.figure()
camera = Camera(fig)

n_particles = 100
box_width = 1000
n_steps = 100
dt = 1
tnow=0

nd=[]

bb=[Cell(0) for i in range(10)]
for j in range(n_steps):
    print(j)
    tnow=dt*j
    
    bb=[b.proliferate() if tnow > b.nextdiv else [b] for b in bb]
    
    bb=sum(bb,[])
    
    


    
    for i in bb:
        kk=i.ssa()
        
        i.xpos += i.xvel*dt
        i.ypos += i.yvel*dt
    
        if abs(i.xpos) > box_width:
          i.xvel = -i.xvel
          i.xpos += i.xvel*dt
        
        
        if abs(i.ypos) > box_width:
          i.yvel = -i.yvel
          i.ypos += i.yvel*dt
    
    x_coord=[b.xpos for b in bb]
    y_coord=[b.ypos for b in bb]
    
    nd.append([x_coord,y_coord])
    
    if j%10 == 0:
        plt.plot(nd[j][0],nd[j][1],'ro')
        camera.snap()
        
animation = camera.animate()
animation.save('celluloid_minimal.gif', writer = 'imagemagick') 
      
       
        # print([i.xpos,i.ypos,i.xvel,i.yvel])
    