
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 21, 2021

@author: Arnab Bandyopadhyay
"""

import numpy as np
import random
import gillespie as gp
from gillespy2.solvers.cpp import SSACSolver
import params
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

from scipy.stats import truncnorm

def get_truncated_normal(mean, sd, low, upp):
    return truncnorm(
        (low - mean) / sd, (upp - mean) / sd, loc=mean, scale=sd)
RN = get_truncated_normal(mean=0.5, sd=0.1, low=0, upp=1)

sigma=params.sigma
box_width = params.box_width


class Cell:
    tnow=0

    def __init__(self, birthtime, gmax, lmt, factor, myx, **kwargs):
        
        
        self.sigma=params.sigma
        self.birthtime=birthtime
        self.gmax=gmax
        self.lmt=lmt
        self.volume=1 
        self.myx=myx
        self.influx=myx

        # self.nextdiv=self.birthtime+np.random.randint(mu-sigma, mu+sigma)#np.random.normal(mu, sigma)
        # self.factor=factor
        # self.generation = 0                                            
        [self.ppd, self.ppd_r, self.ppk, 
        self.ppk_r, self.ppr, self.pg3p, self.pmd,  
        self.pmr, self.pmk, self.pd, 
        self.pr, self.pk, self.pg3pr, self.ppar] = gp.gly_dat[random.randrange(0,199),:]
        self.mu=self.gmax
        self.nextdiv=self.birthtime+self.mu
        self.factor=1 # np.exp(-(self.mu-self.gmax)/50)
        
        if 'xpos' in kwargs: self.xpos = kwargs['xpos']
        else: self.xpos = np.random.randint(box_width) #round(np.random.random(),2)*box_width
        
        if 'ypos' in kwargs: self.ypos = kwargs['ypos']
        else: self.ypos =np.random.randint(box_width) #round(np.random.random(),2)*box_width
        
        if 'xvel' in kwargs: self.xvel = kwargs['xvel']
        else: self.xvel = 2*(np.random.random()-0.5)/2
        
        if 'yvel' in kwargs: self.yvel = kwargs['yvel']
        else: self.yvel =2*(np.random.random()-0.5)/2
    
    
    @classmethod  
    def clone(cls, b):

        if params.inheritance==1:
            return cls(birthtime=b.tnow, gmax=b.gmax, mu=b.mu, myx=b.myx, lmt=b.lmt, volume=1, nextdiv=b.tnow,
                    factor=b.factor, #np.random.choice(params.factors)
                    influx=b.influx,
                    xpos = b.xpos,
                    ypos = b.ypos,
                    xvel=2*(np.random.random()-0.5)/2,
                    yvel=2*(np.random.random()-0.5)/2)
        elif params.inheritance==0:
            nd=np.random.choice(params.myx)
            return cls(birthtime=b.tnow, gmax=b.gmax, mu=b.mu, myx=nd, lmt=b.lmt, volume=1, nextdiv=b.tnow+b.mu,
                    factor=b.factor, #np.random.choice(params.factors)
                    influx=nd,
                    xpos = b.xpos,
                    ypos = b.ypos,
                    xvel=2*(np.random.random()-0.5)/2,
                    yvel=2*(np.random.random()-0.5)/2)
       
    
    def gen(self):
        self.generation += 1

    def proliferate(self, switch, prob, sw_back, reten_prob):
        
        ppd_i=self.ppd
        ppd_r_i=self.ppd_r
        ppk_i=self.ppk
        ppk_r_i=self.ppk_r
        ppr_i=self.ppr
        
        g3p=1*self.pg3p
        md=1*self.pmd
        mr=1*self.pmr
        mk=1*self.pmk
        d=1*self.pd       
        r=1*self.pr
        k=1*self.pk
        g3pr=1*self.pg3pr
        ppar=1*self.ppar
        
        rn=RN.rvs()
        
        g3p1=int(rn*g3p)
        md1=int(rn*md)
        mr1=int(rn*mr)
        mk1=int(rn*mk)
        d1=int(rn*d)
        r1=int(rn*r)
        k1=int(rn*k)
        g3pr1=int(rn*g3pr)
           
        nl=[Cell.clone(self),
            Cell.clone(self)]
        
        if params.inheritance==1:
            if self.influx>=0.4:
                nl[0].influx=self.influx
                if random.uniform(0, 1)<reten_prob:
                    nl[1].influx=random.uniform(0.1,0.2)#self.influx
                else:nl[1].influx=self.influx
            # else:
            #      if ((switch=='yes') and (random.uniform(0, 1)<prob) and (self.influx<=0.2)):
            #          nl[0].influx=random.uniform(0.4, 0.5) # for ngr jumping to 0.4-0.5, for gr 0.5-0.6
            #          nl[1].influx=random.uniform(0.4, 0.5)#self.influx
        
        if ((switch=='yes') and (random.uniform(0, 1)<prob) and (self.influx<=0.2)):
            nl[0].influx=random.uniform(0.4, 0.5) # for ngr jumping to 0.4-0.5, for gr 0.5-0.6
            nl[1].influx=random.uniform(0.4, 0.5)#self.influx
        
        nl[0].ppd, nl[1].ppd = ppd_i,ppd_i
        nl[0].ppd_r,nl[1].ppd_r=ppd_r_i,ppd_r_i
        nl[0].ppk,nl[1].ppk=ppk_i,ppk_i
        nl[0].ppk_r,nl[1].ppk_r=ppk_r_i,ppk_r_i
        nl[0].ppr,nl[1].ppr=ppr_i,ppr_i
        nl[0].pg3p,nl[1].pg3p=g3p1,g3p-g3p1
        nl[0].pmd,nl[1].pmd=md1,md-md1
        nl[0].pmr,nl[1].pmr=mr1, mr-mr1
        nl[0].pmk,nl[1].pmk=mk1, mk-mk1
        nl[0].pd,nl[1].pd=d1, d-d1
        nl[0].pr,nl[1].pr=r1,r-r1
        nl[0].pk,nl[1].pk=k1,k-k1
        nl[0].pg3pr,nl[1].pg3pr=g3pr1,g3pr-g3pr1
        nl[0].ppar, nl[1].ppar = ppar,ppar
        
        nl[0].influx=nl[0].influx+nl[0].pg3p/(params.km+nl[0].pg3p)
        nl[1].influx=nl[1].influx+nl[1].pg3p/(params.km+nl[1].pg3p)
        
        
        nl[0].gmax=np.random.choice(params.gmax)
        nl[0].mu=nl[0].gmax
        nl[0].nextdiv=nl[0].birthtime+nl[0].mu
        
        nl[1].gmax=np.random.choice(params.gmax)
        nl[1].mu=nl[1].gmax
        nl[1].nextdiv=nl[1].birthtime+nl[1].mu
        
       
        
        nl[0].ssa()
        nl[1].ssa()

        return nl
    
    def ssa(self):

        model = gp.Gillespie(self.ppd, self.ppd_r, self.ppk, self.ppk_r,self.ppr, self.pg3p, self.pmd, 
                                      self.pmr, self.pmk, self.pd, self.pr, self.pk, self.pg3pr, self.ppar, self.factor, self.mu, self.influx)
        # results = model.run(TauLeapingSolver)
        results = model.run(NumPySSASolver)
        self.ppd=results['pD'][-1]
        self.ppd_r=results['pD_R'][-1] 
        self.ppk=results['pK'][-1] 
        self.ppk_r=results['pK_R'][-1]
        self.ppr=results['pR'][-1] 
        self.pg3p=results['G3P'][-1]
        self.pmd=results['mD'][-1]  
        self.pmr=results['mR'][-1] 
        self.pmk=results['mK'][-1] 
        self.pd=results['D'][-1] 
        self.pr=results['R'][-1] 
        self.pk=results['K'][-1] 
        self.pg3pr=results['G3PR'][-1] 
        self.ppar=results['paR'][-1] 