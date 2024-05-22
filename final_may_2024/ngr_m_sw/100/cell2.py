
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 21 10:58:43 2021

@author: Arnab
"""

#import matplotlib.pyplot as plt
#from matplotlib.animation import FuncAnimation
#from matplotlib import animation
import numpy as np
import random
#from celluloid import Camera
#import numpy as np
import gillespie as gp
from gillespy2.solvers.cpp import SSACSolver
import scipy.spatial as spatial
from scipy.stats import kde
from numba import jit
import params
import ode_int
#import sde
#import sdeint
from scipy.integrate import odeint,solve_ivp
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

# mu, sigma = 120, 3 
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
        self.influx=myx#random.uniform(0.01,self.myx)
        # self.nextdiv=self.birthtime+np.random.randint(mu-sigma, mu+sigma)#np.random.normal(mu, sigma)
        
        # self.factor=factor
        
            
    
        # self.generation = 0
        
        
        # #### deterministic ode
        # [self.ppd, self.ppd_r, self.pfk, self.pfk_r, self.ppt, 
        # self.ppt_r, self.ppr, self.pg3p, self.pmd, self.pmk, 
        # self.pmr, self.pmf, self.pmt, self.pd, self.pk, 
        # self.pr, self.pf, self.pt, self.pgly, self.pglyk, 
        # self.pg3pr, self.pg3pd] = [1,0,1,0,1,
        # 0,1,np.random.randint(0,lmt),np.random.randint(0,lmt),np.random.randint(0,lmt),
        # np.random.randint(0,lmt),np.random.randint(0,lmt),np.random.randint(0,lmt),np.random.randint(0,lmt),np.random.randint(0,lmt),
        # np.random.randint(0,lmt),np.random.randint(0,lmt),np.random.randint(0,lmt),np.random.randint(0,lmt),np.random.randint(0,lmt),
        # np.random.randint(0,lmt),np.random.randint(0,lmt)]
        
        
        
        ###### simple ode
        
        # [self.ppd, self.ppd_r, self.ppt, 
        # self.ppt_r, self.ppr, self.pg3p, self.pmd,  
        # self.pmr, self.pmt, self.pd, 
        # self.pr, self.pt] = [1,0,1,
        # 0,1,np.random.randint(lmt-5,lmt),np.random.randint(lmt-5,lmt),
        # 2,#np.random.randint(lmt-5,lmt),
        # np.random.randint(lmt-5,lmt),np.random.randint(lmt-5,lmt),
        # 10,#np.random.randint(lmt-5,lmt),
        # np.random.randint(lmt-5,lmt)]
        
        # [self.ppd, self.ppd_r, self.ppk, 
        # self.ppk_r, self.ppr, self.pg3p, self.pmd,  
        # self.pmr, self.pmk, self.pd, 
        # self.pr, self.pk, self.pg3pr, self.ppar] = [1,0,1,
        # 0,1,np.random.randint(0,50),np.random.randint(lmt-5,lmt),
        # 2,#np.random.randint(lmt-5,lmt),
        # np.random.randint(lmt-5,lmt),np.random.randint(lmt-5,lmt),
        # 5,#np.random.randint(lmt-5,lmt),
        # np.random.randint(lmt-5,lmt),np.random.randint(lmt-5,lmt),1]
                                                    
        [self.ppd, self.ppd_r, self.ppk, 
        self.ppk_r, self.ppr, self.pg3p, self.pmd,  
        self.pmr, self.pmk, self.pd, 
        self.pr, self.pk, self.pg3pr, self.ppar] = gp.gly_dat[random.randrange(0,199),:]
        
        # self.pd=np.random.choice([1,2,3,4,5,1,3,5,60,70,80])
                                                    
                                                    
        # [self.ppd, self.ppd_r, self.ppk, 
        # self.ppk_r, self.ppr, self.pg3p, self.pmd,  
        # self.pmr, self.pmk, self.pd, 
        # self.pr, self.pk, self.pg3pr, self.ppar] = [1,0,1,0,1,10,1,1,1,1,1,1,1,1]
                                                    
                                                    
        # self.mu=90/(1+self.pd/(20+self.pd))
        # self.mu=30*(1+2*(self.pd*self.pd)/(500+self.pd*self.pd))
        # self.mu=40*(1+2*(self.pr*self.pr)/(300+self.pr*self.pr))
        # if self.pr>=20:
        #     # self.mu=params.gmax*(np.exp(0.7*self.pr*self.pr/(100+self.pr*self.pr)))
        #     self.mu=random.randrange(50,100)
        # else:
        #     self.mu=random.randrange(30,40)
        
        # self.mu=np.random.choice(params.gmax)*(np.exp(0.9*self.pr*self.pr/(100+self.pr*self.pr)))#0.7-1.0
        # self.mu=self.gmax*(np.exp(0.5*self.pr*self.pr/(100+self.pr*self.pr)))#0.7-1.0
        # self.mu=self.gmax*(1+1.2*self.pr*self.pr/(150*150+self.pr*self.pr))#0.7-1.0
#        self.mu=self.gmax*(1+3*self.pr*self.pr/(100+self.pr*self.pr))
        self.mu=self.gmax
        
        # if self.pr>=20:
        #     self.mu=self.gmax*(1+1.5*self.pr*self.pr/(100+self.pr*self.pr))
        # else:
        #     self.mu=self.gmax
        
#        if self.mu<30:
#            self.mu=30
        self.nextdiv=self.birthtime+self.mu
        self.factor=1 # np.exp(-(self.mu-self.gmax)/50)
        
        # self.mu=30*(1+2*self.pg3p/(20+self.pg3p))
                                         
        # [self.ppd, self.ppd_r, self.ppk, 
        # self.ppk_r, self.ppr, self.pg3p, self.pmd,  
        # self.pmr, self.pmk, self.pd, 
        # self.pr, self.pk, self.pg3pr] = [1,0,1,
        # 0,1,np.random.randint(0,lmt),np.random.randint(0,lmt),
        # np.random.randint(0,lmt),
        # np.random.randint(0,lmt),np.random.randint(0,lmt),
        # np.random.randint(0,lmt),
        # np.random.randint(0,lmt),np.random.randint(0,lmt)]
        
        if 'xpos' in kwargs: self.xpos = kwargs['xpos']
        else: self.xpos = np.random.randint(box_width) #round(np.random.random(),2)*box_width
        
        if 'ypos' in kwargs: self.ypos = kwargs['ypos']
        else: self.ypos =np.random.randint(box_width) #round(np.random.random(),2)*box_width
        
        if 'xvel' in kwargs: self.xvel = kwargs['xvel']
        else: self.xvel = 2*(np.random.random()-0.5)/2
        
        if 'yvel' in kwargs: self.yvel = kwargs['yvel']
        else: self.yvel =2*(np.random.random()-0.5)/2
        
        
        
            
    # ##### deterministic ode    
    # @classmethod  
    # def clone(cls, b,g3p,md,mk,mr,mf,mt,d,k,r,f,t,gly,glyk,g3pr,g3pd):
    #     return cls(birthtime=b.tnow, mu=b.mu, lmt=b.lmt, volume=1, nextdiv=b.tnow+np.random.randint(b.mu-sigma,b.mu+sigma),
    #                ppd = 1, ppd_r=0, pfk = 1, pfk_r=0, ppt=1, ppt_r=0, ppr=0,
    #                pg3p=g3p, pmd=md,pmk=mk,pmr=mr,pmf=mf, pmt=mt,pd=d,pk=k,pr=r, pf=f, pt=t, pgly=gly, 
    #                pglyk=glyk, pg3pr=g3pr, pg3pd=g3pd,
    #                xpos = b.xpos,
    #                ypos = b.ypos,
    #                xvel=2*(np.random.random()-0.5)/2,
    #                yvel=2*(np.random.random()-0.5)/2)
    
     ##### simple ode    
    # @classmethod  
    # def clone(cls, b,ppd_i,ppd_r_i,ppt_i,ppt_r_i,ppr_i,g3p,md,mr,mt,d,r,t):
    #     return cls(birthtime=b.tnow, mu=b.mu, lmt=b.lmt, volume=1, nextdiv=b.tnow+np.random.randint(b.mu-sigma,b.mu+sigma),
    #                ppd = ppd_i, ppd_r=ppd_r_i, ppt=ppt_i, ppt_r=ppt_r_i, ppr=ppr_i,
    #                pg3p=g3p, pmd=md, pmr=mr, pmt=mt,pd=d,pr=r, pt=t,
    #                # factor=b.factor,
    #                factor=np.random.choice(params.factors),
    #                xpos = b.xpos,
    #                ypos = b.ypos,
    #                xvel=2*(np.random.random()-0.5)/2,
    #                yvel=2*(np.random.random()-0.5)/2)
    
    
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
            # ndiv=np.random.randint(b.mu-b.sigma, b.mu+b.sigma)
        # if ndiv<30:
        #     ndiv=30
        # ndiv=30*(1+2*b.pd/(20+b.pd))
        # ndiv=90/(1+b.pd/(20+b.pd))
        # ndiv=30*(1+2*b.pg3p/(20+b.pg3p))
        # ndiv=30*(1+2*(b.pd*b.pd)/(500+b.pd*b.pd))
            return cls(birthtime=b.tnow, gmax=b.gmax, mu=b.mu, myx=nd, lmt=b.lmt, volume=1, nextdiv=b.tnow+b.mu,
                    factor=b.factor, #np.random.choice(params.factors)
                    influx=nd,
                    xpos = b.xpos,
                    ypos = b.ypos,
                    xvel=2*(np.random.random()-0.5)/2,
                    yvel=2*(np.random.random()-0.5)/2)
       
        
    
    
    
    
    # def clone(self,ppd_i,ppd_r_i,ppt_i,ppt_r_i,ppr_i,g3p,md,mr,mt,d,r,t):
    #     birthtime=self.tnow
    #     mu=self.mu
    #     lmt=self.lmt
    #     volume=1
    #     nextdiv=self.tnow+np.random.randint(self.mu-sigma,self.mu+sigma)
    #     ppd = ppd_i
    #     ppd_r=ppd_r_i
    #     ppt=ppt_i
    #     ppt_r=ppt_r_i
    #     ppr=ppr_i
    #     pg3p=g3p
    #     pmd=md
    #     pmr=mr
    #     pmt=mt
    #     pd=d
    #     pr=r 
    #     pt=t
    #     # factor=b.factor,
    #     factor=np.random.choice(params.factors)
    #     xpos = self.xpos
    #     ypos = self.ypos
    #     xvel=2*(np.random.random()-0.5)/2
    #     yvel=2*(np.random.random()-0.5)/2
        
    def gen(self):
        self.generation += 1

            
        ####### simple ode section
    def proliferate(self, switch, prob, sw_back, reten_prob):
        # newborn_mu1=np.random.randint(self.mu-self.sigma,self.mu+self.sigma)
        # newborn_mu2=np.random.randint(self.mu-self.sigma,self.mu+self.sigma)
        # if newborn_mu1<30:
        #     newborn_mu1=30
        # if newborn_mu2<30:
        #     newborn_mu2=30
        
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
        
        # if np.isnan(g3p)==True:
        #     g3p=0
            
        # if np.isnan(md)==True:
        #     md=0
        # if np.isnan(mr)==True:
        #     mr=0
        # if np.isnan(mk)==True:
        #     mk=0
        # if np.isnan(d)==True:
        #     d=0
        # if np.isnan(r)==True:
        #     r=0
        # if np.isnan(k)==True:
        #     k=0
       
        
        # g3p1=int(0.7*g3p)#int(np.random.random()*g3p)
        # md1=int(0.7*md)#int(np.random.random()*md)
        # mk1=int(0.7*mk)#int(np.random.random()*mk)
        # mr1=int(0.7*mr)#int(np.random.random()*mr)
        # d1=int(0.7*d)#int(np.random.random()*d)
        # k1=int(0.7*k)#int(np.random.random()*k)
        # r1=int(0.7*r)#int(np.random.random()*r)
        # gly1=int(0.7*gly)#int(np.random.random()*gly)
        # g3pr1=int(0.7*g3pr)#int(np.random.random()*g3pr)
        
#        g3p1=int(np.random.random()*g3p)
#        md1=int(np.random.random()*md)
#        mr1=int(np.random.random()*mr)
#        mk1=int(np.random.random()*mk)
#        d1=int(np.random.random()*d)
#        r1=int(np.random.random()*r)
#        k1=int(np.random.random()*k)
#        g3pr1=int(np.random.random()*g3pr)
        
        rn=RN.rvs()
        
        g3p1=int(rn*g3p)
        md1=int(rn*md)
        mr1=int(rn*mr)
        mk1=int(rn*mk)
        d1=int(rn*d)
        r1=int(rn*r)
        k1=int(rn*k)
        g3pr1=int(rn*g3pr)
        
        
        
        # g3p1=int(0.9*g3p)
        # md1=int(0.9*md)
        # mr1=int(0.9*mr)
        # mk1=int(0.9*mk)
        # d1=int(0.9*d)
        # r1=int(0.9*r)
        # k1=int(0.9*k)
        # g3pr1=int(0.9*g3pr)
       
#        print('before switching')
#        print(self.influx)
        
        # nl=[Cell.clone(self,ppd_i,ppd_r_i,ppt_i,ppt_r_i,ppr_i,g3p1, md1, mr1, mt1, d1, r1, t1),
        #     Cell.clone(self,ppd_i,ppd_r_i,ppt_i,ppt_r_i,ppr_i,g3p-g3p1,md-md1,mr-mr1, mt-mt1, 
        #                 d-d1,r-r1, t-t1)]
#        
#        if switch=='yes':
##            print('switching')
#            if ((random.uniform(0, 1)<prob) and (self.influx<=0.3)):
##                print('switching')
#                
#                self.influx=random.uniform(0.5, 0.6)
#                print(self.influx)
          
        nl=[Cell.clone(self),
            Cell.clone(self)]
        
        if params.inheritance==1:
            if self.influx>=0.4:
                nl[0].influx=self.influx
                if random.uniform(0, 1)<reten_prob:
                    nl[1].influx=random.uniform(0.1,0.2)#self.influx
                else:nl[1].influx=self.influx
            else:
                 if ((switch=='yes') and (random.uniform(0, 1)<prob) and (self.influx<=0.2)):
                     nl[0].influx=random.uniform(0.4, 0.5) # for ngr jumping to 0.4-0.5, for gr 0.5-0.6
                     nl[1].influx=random.uniform(0.4, 0.5)#self.influx
        
        # if ((switch=='yes') and (random.uniform(0, 1)<prob) and (self.influx<=0.2)):
        #     nl[0].influx=random.uniform(0.4, 0.5) # for ngr jumping to 0.4-0.5, for gr 0.5-0.6
        #     nl[1].influx=random.uniform(0.4, 0.5)#self.influx






#        
#        nl[0].influx=self.influx
#        nl[1].influx=random.uniform(0.2, 0.3)#self.influx
#        
#        if switch=='yes':
#            if ((random.uniform(0, 1)<sw_back) and (nl[1].influx>0.3)):
#                nl[1].influx=random.uniform(0.2, 0.3)
#            if ((random.uniform(0, 1)<sw_back) and (nl[0].influx>0.3)):
#                nl[0].influx=random.uniform(0.2, 0.3)
                    
                
           
               
           
               
        
       
        
           
           
        
               
              
        
#        print('new born influx')
#        print([nl[0].influx, nl[1].influx])
        # nl=[Cell.clone(self,g3p, md, mk, mr, d, k, r, gly, g3pr),
        #     Cell.clone(self,g3p, md, mk, mr, d, k, r, gly, g3pr)]
        
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
        
        nl[0].influx=nl[0].influx+nl[0].pg3p/(100+nl[0].pg3p)
        nl[1].influx=nl[1].influx+nl[1].pg3p/(100+nl[1].pg3p)
        
        # [nl[0].ppd, nl[0].ppd_r, nl[0].ppk, 
        # nl[0].ppk_r, nl[0].ppr, nl[0].pg3p, nl[0].pmd,  
        # nl[0].pmr, nl[0].pmk, nl[0].pd, 
        # nl[0].pr, nl[0].pk, nl[0].pg3pr, nl[0].ppar] = [1,0,1,0,1,10,1,1,1,1,1,1,1,1]
        
        # [nl[1].ppd, nl[1].ppd_r, nl[1].ppk, 
        # nl[1].ppk_r, nl[1].ppr, nl[1].pg3p, nl[1].pmd,  
        # nl[1].pmr, nl[1].pmk, nl[1].pd, 
        # nl[1].pr, nl[1].pk, nl[1].pg3pr, nl[1].ppar] = [1,0,1,0,1,10,1,1,1,1,1,1,1,1]
        
#        if params.inheritance==0:
            # if nl[0].pr>=20:
            #     nl[0].mu=params.gmax*(np.exp(nl[0].pr*nl[0].pr/(200+nl[0].pr*nl[0].pr)))
            # else:
            #     nl[0].mu=random.randrange(30,60)
                
            # if nl[1].pr>=20:
            #     nl[1].mu=params.gmax*(np.exp(nl[1].pr*nl[1].pr/(200+nl[1].pr*nl[1].pr)))
            # else:
            #     nl[1].mu=random.randrange(30,60)
            
            
            # ndiv1=np.random.choice(params.gmax)*(np.exp(0.7*nl[0].pr*nl[0].pr)/(1000+nl[0].pr*nl[0].pr))
            # ndiv2=np.random.choice(params.gmax)*(np.exp(0.7*nl[1].pr*nl[1].pr)/(1000+nl[1].pr*nl[1].pr))
            
#            ndiv1=nl[0].gmax*(1+1.2*nl[0].pr*nl[0].pr/(150*150+nl[0].pr*nl[0].pr))#0.7-1.0
#            ndiv2=nl[1].gmax*(1+1.2*nl[1].pr*nl[1].pr/(150*150+nl[1].pr*nl[1].pr))#0.7-1.0
#            if ndiv1<=30:
#                ndiv1=30
#            
#            if ndiv2<=30:
#                ndiv2=30
#            
#            nl[0].mu, nl[1].mu = ndiv1, ndiv2
#            nl[0].nextdiv, nl[1].nextdiv = nl[0].tnow+ndiv1, nl[1].tnow+ndiv2
#            
#            nl[0].factor=np.exp(-(nl[0].mu-30)/50)
#            nl[1].factor=np.exp(-(nl[1].mu-30)/50)
            # if ndiv1>=40:
            #     nl[0].factor=np.exp(-(ndiv1-30)/50)
            # else:
            #     nl[0].factor=1
            
            # if ndiv2>=40:
            #     nl[1].factor=np.exp(-(ndiv2-30)/50)
            # else:
            #     nl[1].factor=1
    
        
        # nl[0].gmax=20+30*np.exp(-0.1*nl[0].pd)
        # nl[0].mu=nl[0].gmax
        # nl[0].nextdiv=nl[0].birthtime+nl[0].mu
        
        # nl[1].gmax=20+30*np.exp(-0.1*nl[1].pd)
        # nl[1].mu=nl[1].gmax
        # nl[1].nextdiv=nl[1].birthtime+nl[1].mu
        
        nl[0].gmax=np.random.choice(params.gmax)
        nl[0].mu=nl[0].gmax
        nl[0].nextdiv=nl[0].birthtime+nl[0].mu
        
        nl[1].gmax=np.random.choice(params.gmax)
        nl[1].mu=nl[1].gmax
        nl[1].nextdiv=nl[1].birthtime+nl[1].mu
        
       
        
        nl[0].ssa()
        nl[1].ssa()
        
        
        
       
        
#        if nl[0].pd>15:
#            nl[0].gmax=20
#            nl[0].mu=20
#            nl[0].nextdiv=nl[0].birthtime+nl[0].mu
#        else:
#            nl[0].gmax=np.random.choice(params.gmax)
#            nl[0].mu=nl[0].gmax
#            nl[0].nextdiv=nl[0].birthtime+nl[0].mu
#                
#        if nl[1].pd>15:
#            nl[1].gmax=20
#            nl[1].mu=20
#            nl[1].nextdiv=nl[1].birthtime+nl[1].mu
#        else:
#            nl[1].gmax=np.random.choice(params.gmax)
#            nl[1].mu=nl[1].gmax
#            nl[1].nextdiv=nl[1].birthtime+nl[1].mu

        return nl
        
        
        
        
        
        
     #### Deterministic ode section       

    # def proliferate(self):

    #     g3p=self.pg3p
    #     md=self.pmd
    #     mk=self.pmk
    #     mr=self.pmr
    #     mf=self.pmf
    #     mt=self.pmt
    #     d=self.pd
    #     k=self.pk
    #     r=self.pr
    #     f=self.pf
    #     t=self.pt
    #     gly=self.pgly
    #     glyk=self.pglyk
    #     g3pr=self.pg3pr
    #     g3pd=self.pg3pd
        
    #     # g3p1=int(0.7*g3p)#int(np.random.random()*g3p)
    #     # md1=int(0.7*md)#int(np.random.random()*md)
    #     # mk1=int(0.7*mk)#int(np.random.random()*mk)
    #     # mr1=int(0.7*mr)#int(np.random.random()*mr)
    #     # d1=int(0.7*d)#int(np.random.random()*d)
    #     # k1=int(0.7*k)#int(np.random.random()*k)
    #     # r1=int(0.7*r)#int(np.random.random()*r)
    #     # gly1=int(0.7*gly)#int(np.random.random()*gly)
    #     # g3pr1=int(0.7*g3pr)#int(np.random.random()*g3pr)
        
    #     g3p1=int(np.random.random()*g3p)
    #     md1=int(np.random.random()*md)
    #     mk1=int(np.random.random()*mk)
    #     mr1=int(np.random.random()*mr)
    #     mf1=int(np.random.random()*mf)
    #     mt1=int(np.random.random()*mt)
    #     d1=int(np.random.random()*d)
    #     k1=int(np.random.random()*k)
    #     r1=int(np.random.random()*r)
    #     f1=int(np.random.random()*f)
    #     t1=int(np.random.random()*t)
    #     gly1=int(np.random.random()*gly)
    #     glyk1=int(np.random.random()*glyk)
    #     g3pr1=int(np.random.random()*g3pr)
    #     g3pd1=int(np.random.random()*g3pd)
        
        
    #     nl=[Cell.clone(self,g3p1,md1, mk1, mr1,mf1, mt1, d1, k1, r1, f1, t1, gly1, glyk1, g3pr1, g3pd1),
    #         Cell.clone(self,g3p-g3p1,md-md1,mk-mk1,mr-mr1,mf-mf1, mt-mt1, 
    #                    d-d1,k-k1,r-r1, f-f1, t-t1, gly-gly1, glyk-glyk1, g3pr-g3pr1, g3pd-g3pd1)]
        
    #     # nl=[Cell.clone(self,g3p, md, mk, mr, d, k, r, gly, g3pr),
    #     #     Cell.clone(self,g3p, md, mk, mr, d, k, r, gly, g3pr)]

    #     return nl
    
    # def ssa(self):
    #     ita=self.pd/(70+self.pd)
    #     model = gp.Gillespie(self.pinpd,self.pinpk,self.pactpd,
    #                                self.pactpk,self.ppr, self.pg3p,self.pmd, self.pmk, self.pmr,
    #                                self.pd, self.pk, self.pr, self.pgly, self.pg3pr, self.ppar, ita)
    #     results = model.run(TauLeapingSolver)
    #     # results = model.run(NumPySSASolver)
    #     self.pinpd=results['inPD_R'][-1]
    #     self.pinpk=results['inPK_R'][-1]
    #     self.pactpd=results['actPD'][-1]
    #     self.pactpk=results['actPK'][-1]
    #     self.ppr=results['PR'][-1]
    #     self.pg3p=results['G3P'][-1]
    #     self.pmd=results['mD'][-1]
    #     self.pmk=results['mK'][-1]
    #     self.pmr=results['mR'][-1]
    #     self.pd=results['D'][-1]
    #     self.pk=results['K'][-1]
    #     self.pr=results['R'][-1]
    #     self.pgly=results['Gly'][-1]
    #     self.pg3pr=results['G3PR'][-1]
    #     self.ppar=results['PAR'][-1]
    
    def ssa(self):
        # ita=self.pd/(70+self.pd)
        # model = gp.Gillespie_bistable(self.ppd, self.ppd_r, self.pfk, self.pfk_r, self.ppt, self.ppt_r,
        #           self.ppr, self.pg3p, self.pmd, self.pmk, self.pmr, self.pmf, self.pmt, self.pd, self.pk, 
        #           self.pr, self.pf, self.pt, 
        #           self.pgly, self.pglyk, self.pg3pr, self.pg3pd, ita)
       
        # print('running ssa')
        model = gp.Gillespie(self.ppd, self.ppd_r, self.ppk, self.ppk_r,self.ppr, self.pg3p, self.pmd, 
                                      self.pmr, self.pmk, self.pd, self.pr, self.pk, self.pg3pr, self.ppar, self.factor, self.mu, self.influx)
        # results = model.run(TauLeapingSolver)
        results = model.run(NumPySSASolver)
        self.ppd=results['pD'][-1]
        self.ppd_r=results['pD_R'][-1] 
        # self.pfk=results['pFK'][-1]
        # self.pfk_r=results['pFK_R'][-1] 
        self.ppk=results['pK'][-1] 
        self.ppk_r=results['pK_R'][-1]
        self.ppr=results['pR'][-1] 
        self.pg3p=results['G3P'][-1]
        self.pmd=results['mD'][-1] 
        # self.pmk=results['mK'][-1] 
        self.pmr=results['mR'][-1] 
        # self.pmf=results['mF'][-1] 
        self.pmk=results['mK'][-1] 
        self.pd=results['D'][-1] 
        # self.pk=results['K'][-1] 
        self.pr=results['R'][-1] 
        # self.pf=results['F'][-1] 
        self.pk=results['K'][-1] 
        self.pg3pr=results['G3PR'][-1] 
        self.ppar=results['paR'][-1] 
        # self.pgly=results['Gly'][-1] 
        # self.pglyk=results['GlyK'][-1] 
        # self.pg3pr=results['G3PR'][-1] 
        # self.pg3pd=results['G3PD'][-1] 
        
    def deterministic_ode(self):
            
        # t = np.linspace(0,5,100)
        z0=[self.ppd_r, self.pfk_r, self.ppt_r,
             self.pmd, self.pmf,  self.pmk, self.pmt, self.pmr, self.pd, self.pf, self.pk, self.pt, self.pr, 
             self.pgly, self.pg3p, self.pg3pr, self.pglyk,  self.pg3pd]
        results = solve_ivp(ode_int.model,(0,2),z0, method='LSODA',t_eval=np.linspace(0,2,100).tolist())
        self.ppd=1-results.y[0,-1]
        self.ppd_r=results.y[0,-1] 
        self.pfk=1-results.y[1,-1]
        self.pfk_r=results.y[1,-1] 
        self.ppt=1-results.y[2,-1]
        self.ppt_r=results.y[2,-1]
        self.ppr=1
        self.pg3p=results.y[3,-1]
        self.pmd=results.y[4,-1] 
        self.pmk=results.y[5,-1] 
        self.pmr=results.y[6,-1] 
        self.pmf=results.y[7,-1] 
        self.pmt=results.y[8,-1] 
        self.pd=results.y[9,-1] 
        self.pk=results.y[10,-1] 
        self.pr=results.y[11,-1] 
        self.pf=results.y[12,-1] 
        self.pt=results.y[13,-1] 
        self.pgly=results.y[14,-1] 
        self.pglyk=results.y[15,-1] 
        self.pg3pr=results.y[16,-1] 
        self.pg3pd=results.y[17,-1] 
        
    def simple_ode(self):
        
        # t = np.linspace(0,5,100)
        z0=[self.ppd_r, self.ppt_r,
             self.pmd, self.pmt, self.pmr, self.pd, self.pt, self.pr, 
             self.pg3p]
        results = solve_ivp(ode_int.simple_model,(0,2),z0, method='LSODA',t_eval=np.linspace(0,2,20).tolist())
        self.ppd=1-results.y[0,-1]
        self.ppd_r=results.y[0,-1] 
        self.ppt=1-results.y[1,-1]
        self.ppt_r=results.y[1,-1]
        self.ppr=1
        self.pmd=results.y[2,-1] 
        self.pmt=results.y[3,-1]
        self.pmr=results.y[4,-1] 
        self.pd=results.y[5,-1] 
        self.pt=results.y[6,-1] 
        self.pr=results.y[7,-1] 
        self.pg3p=results.y[8,-1]
        
        
    def simple_sde(self):
        
        # t = np.linspace(0,5,100)
        z0=[self.ppd_r, self.ppt_r,
             self.pmd, self.pmt, self.pmr, self.pd, self.pt, self.pr, 
             self.pg3p]
        tspan = np.linspace(0, 2., 1000)
        results = sdeint.itoEuler(sde.simple_model, sde.GG, z0, tspan)
        self.ppd=1-results[-1,0]
        self.ppd_r=results[-1,0] 
        self.ppt=1-results[-1,1]
        self.ppt_r=results[-1,1]
        self.ppr=1
        self.pmd=results[-1,2] 
        self.pmt=results[-1,3]
        self.pmr=results[-1,4] 
        self.pd=results[-1,5] 
        self.pt=results[-1,6] 
        self.pr=results[-1,7] 
        self.pg3p=results[-1,8]
            
        
        

x_vel, y_vel= 0.5, 0.5



# def take_step(b):
    
#     xvel=b[3]  
#     b[0] += xvel*dt
#     b[1] += b[4]*dt
    
#     if abs(b[0]) > box_width:
#       b[3] = -b[3]
#       b[0] += b[3]*dt
    
#     if abs(b[1]) > box_width:
#       b[4] = -b[4]
#       b[1] += b[4]*dt
    
#     return b

# def update(list1,list2,molecule):
#     if molecule=='Gly':
#         stepsize=5
#     elif molecule=='G3P':
#         stepsize=2
        
#     for i in range(len(list1)):
#         if list1[i] == "RIGHT":
#             if list2[i,0]<box_width:
#                 list2[i,0]+=stepsize
#         elif list1[i] == "LEFT":
#             if list2[i,0]>0:
#                 list2[i,0]-=stepsize
#         elif list1[i] == "UP":
#             if list2[i,1]<box_width:
#                 list2[i,1]+=stepsize
#         elif list1[i] == "DOWN":
#             if list2[i,1] > 0:
#                 list2[i,1]-=stepsize
#     return list2

@jit(nopython=True)
def update_gly(list1,list2):
    if list1 == 1:
        if list2[0]<box_width:
            list2[0]+=4
    elif list1 == -1:
        if list2[0]>0:
            list2[0]-=4
    elif list1 == 10:
        if list2[1]<box_width:
            list2[1]+=4
    elif list1 == -10:
        if list2[1] > 0:
            list2[1]-=4
    return list2

@jit(nopython=True)
def update_g3p(list1,list2):
    if list1 == 1:
        if list2[0]<box_width:
            list2[0]+=2
    elif list1 == -1:
        if list2[0]>0:
            list2[0]-=2
    elif list1 == 10:
        if list2[1]<box_width:
            list2[1]+=2
    elif list1 == -10:
        if list2[1] > 0:
            list2[1]-=2
    return list2


# def update2_gly(list1,list2):
#     if list1 == "RIGHT":
#         if list2[0]<box_width:
#             list2[0]+=5
#     elif list1 == "LEFT":
#         if list2[0]>0:
#             list2[0]-=5
#     elif list1 == "UP":
#         if list2[1]<box_width:
#             list2[1]+=5
#     elif list1 == "DOWN":
#         if list2[1] > 0:
#             list2[1]-=5
#     return list2

# def update2_g3p(list1,list2):
#     if list1 == "RIGHT":
#         if list2[0]<box_width:
#             list2[0]+=2
#     elif list1 == "LEFT":
#         if list2[0]>0:
#             list2[0]-=2
#     elif list1 == "UP":
#         if list2[1]<box_width:
#             list2[1]+=2
#     elif list1 == "DOWN":
#         if list2[1] > 0:
#             list2[1]-=2
#     return list2


@jit(nopython=True)                
def osmotic(inside,cell_volume,outside,radius):
    """ 
    osmotic pressure= i*c*R*T; i, R & T is constant for inside and outside
    cell volume=2.4 um2
    
    """
    ci=inside/(2.4*cell_volume)
    co=outside/(np.pi*radius*radius)
    ch=(inside+outside)/(2.4*cell_volume+np.pi*radius*radius)
    
    return round(ch*2.4*cell_volume)
    


# def uptake(a,b):

#     for i in range(len(b)):
        
#         if g3px_xy[b[i],2]>2:
#             a.pg3p += 2
#             g3px_xy[b[i],2] -= 2
#             break
        


# # fig = plt.figure()
# # camera = Camera(fig)






# nd=[]
# g3px_xy=np.array([[0,0,0,0.1,0.1]])


# bb=[Cell(0) for i in range(5)]
# for j in range(n_steps):
#     print(j)
#     tnow=dt*j
    
#     bb=[b.proliferate() if tnow > b.nextdiv else [b] for b in bb]
    
#     bb=sum(bb,[])
    
    


    
#     for i in bb:
#         kk=i.ssa()
        
        
#         g3px_xy=np.append(g3px_xy, [[i.xpos,i.ypos, kk,
#                                       2*(np.random.random()-0.5),
#                                       2*(np.random.random()-0.5)]],axis=0)
        
        
#         # print([i.xpos,i.ypos])
#         # print([[x,y] for x,y in zip(g3px_x,g3px_y)])
        
#         i.xpos += i.xvel*dt
#         i.ypos += i.yvel*dt
    
#         if abs(i.xpos) > box_width:
#           i.xvel = -i.xvel
#           i.xpos += i.xvel*dt
        
        
#         if abs(i.ypos) > box_width:
#           i.yvel = -i.yvel
#           i.ypos += i.yvel*dt
    
#     if j==0:
#         g3px_xy=np.delete(g3px_xy,0,axis=0)
    
#     # g3px_xy=np.array(g3px_xy)
    
#     g3px_xy=np.array(list(map(take_step,g3px_xy)))
    
#     point_tree = spatial.cKDTree(g3px_xy[:,[0,1]])
    
#     clpt=[point_tree.query_ball_point([b.xpos,b.ypos], 1) for b in bb]
    
#     list(map(uptake,bb,clpt))
    
#     x_coord=[b.xpos for b in bb]
#     y_coord=[b.ypos for b in bb]
    
#     nd.append([x_coord,y_coord])
    
#     # if j%10 == 0:
#     #     plt.plot(nd[j][0],nd[j][1],'ro')
#     #     camera.snap()
        
# #     if j%10 == 0:
# #         xcord=[b[0] for b in g3px_xy]
# #         ycord=[b[1] for b in g3px_xy]
# #         plt.plot(nd[j][0],nd[j][1],'ro')
# #         plt.plot(xcord, ycord,'b.')
# #         camera.snap()
        
# # animation = camera.animate()
# # animation.save('celluloid_minimal_2.gif', writer = 'imagemagick') 
      
       
#         # print([i.xpos,i.ypos,i.xvel,i.yvel])

# # Creating histogram

# hist=[b.pd for b in bb]
# density = kde.gaussian_kde(hist)
# x = np.linspace(min(hist),max(hist),100)
# y=density(x)

# plt.plot(x, y)
# plt.title("Density Plot of the data")
# plt.show()
# print(np.mean(hist))
