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
import scipy.spatial as spatial
from scipy.stats import kde
from numba import jit
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



# mu, sigma = 120, 3 
sigma=params.sigma
box_width = params.box_width


class Cell:
    tnow=0
    
    def __init__(self, birthtime, mu, **kwargs):
        
        
        
        self.birthtime=birthtime
        self.mu=mu
        self.volume=1 
        self.nextdiv=self.birthtime+np.random.randint(mu-sigma, mu+sigma)#np.random.normal(mu, sigma)
        
        if 'pinpd' in kwargs: self.pinpd = kwargs['pinpd']
        else: self.pinpd = 1
        
        if 'pinpk' in kwargs: self.pinpk = kwargs['pinpk']
        else: self.pinpk = 1
        
        if 'pactpd' in kwargs: self.pactpd = kwargs['pactpd']
        else: self.pactpd = 0
        
        if 'pactpk' in kwargs: self.pactpk = kwargs['pactpk']
        else: self.pactpk = 0
        
        if 'ppr' in kwargs: self.ppr = kwargs['ppr']
        else: self.ppr = 1
        
        if 'pg3p' in kwargs: self.pg3p = kwargs['pg3p']
        else: self.pg3p = np.random.randint(0,5)
        
        if 'pmd' in kwargs: self.pmd = kwargs['pmd']
        else: self.pmd = np.random.randint(0,5)
        
        if 'pmk' in kwargs: self.pmk = kwargs['pmk']
        else: self.pmk = np.random.randint(0,5)
        
        if 'pmr' in kwargs: self.pmr = kwargs['pmr']
        else: self.pmr = np.random.randint(0,5)
        
        if 'pd' in kwargs: self.pd = kwargs['pd']
        else: self.pd = np.random.randint(0,5)
        
        if 'pk' in kwargs: self.pk = kwargs['pk']
        else: self.pk = np.random.randint(0,5)
        
        if 'pr' in kwargs: self.pr = kwargs['pr']
        else: self.pr = np.random.randint(0,5)
        
        if 'pgly' in kwargs: self.pgly = kwargs['pgly']
        else: self.pgly = np.random.randint(0,5)
        
        
        # self.generation = 0
        
        if 'xpos' in kwargs: self.xpos = kwargs['xpos']
        else: self.xpos = np.random.randint(box_width) #round(np.random.random(),2)*box_width
        
        if 'ypos' in kwargs: self.ypos = kwargs['ypos']
        else: self.ypos =np.random.randint(box_width) #round(np.random.random(),2)*box_width
        
        if 'xvel' in kwargs: self.xvel = kwargs['xvel']
        else: self.xvel = 2*(np.random.random()-0.5)/2
        
        if 'yvel' in kwargs: self.yvel = kwargs['yvel']
        else: self.yvel =2*(np.random.random()-0.5)/2
        
        
        
            
        
    @classmethod  
    def clone(cls, b,g3p,md,mk,mr,d,k,r,gly):
        return cls(birthtime=b.tnow, mu=b.mu, volume=1, nextdiv=b.tnow+np.random.randint(b.mu-sigma,b.mu+sigma),
                   pinpd = 1, pinpk = 1, pactpd = 0, pactpk = 0, ppr=1,
                   pg3p=g3p, pmd=md,pmk=mk,pmr=mr,pd=d,pk=k,pr=r, pgly=gly,
                   xpos = b.xpos,
                   ypos = b.ypos,
                   xvel=2*(np.random.random()-0.5)/2,
                   yvel=2*(np.random.random()-0.5)/2)
    def gen(self):
        self.generation += 1

    def proliferate(self):

        g3p=self.pg3p
        md=self.pmd
        mk=self.pmk
        mr=self.pmr
        d=self.pd
        k=self.pk
        r=self.pr
        gly=self.pgly
        
        g3p1=int(np.random.random()*g3p)
        md1=int(np.random.random()*md)
        mk1=int(np.random.random()*mk)
        mr1=int(np.random.random()*mr)
        d1=int(np.random.random()*d)
        k1=int(np.random.random()*k)
        r1=int(np.random.random()*r)
        gly1=int(np.random.random()*gly)
        
        
        nl=[Cell.clone(self,g3p1,md1, mk1, mr1, d1, k1, r1, gly1),
            Cell.clone(self,g3p-g3p1,md-md1,mk-mk1,mr-mr1,d-d1,k-k1,r-r1, gly-gly1)]

        return nl
    
    def ssa(self):
        model = gp.MichaelisMenten(self.pinpd,self.pinpk,self.pactpd,
                                   self.pactpk,self.ppr, self.pg3p,self.pmd, self.pmk, self.pmr,
                                   self.pd, self.pk, self.pr, self.pgly)
        results = model.run(TauLeapingSolver)
        self.pinpd=results['inPD_R'][-1]
        self.pinpk=results['inPK_R'][-1]
        self.pactpd=results['actPD'][-1]
        self.pactpk=results['actPK'][-1]
        self.ppr=results['PR'][-1]
        self.pg3p=results['G3P'][-1]
        self.pmd=results['mD'][-1]
        self.pmk=results['mK'][-1]
        self.pmr=results['mR'][-1]
        self.pd=results['D'][-1]
        self.pk=results['K'][-1]
        self.pr=results['R'][-1]
        self.pgly=results['Gly'][-1]
        
        
        # print([results['inPD'][0],results['inPK'][0],
        #        results['actPD'][0],results['actPK'][0],
        #        results['mD'][0],results['mK'][0],
        #        results['D'][0],results['K'][0],
        #        results['G3P'][0],results['G3PX'][0]])
        
        # return results['G3PX'][-1]
        

x_vel, y_vel= 0.5, 0.5


def take_step(b):
    
    xvel=b[3]  
    b[0] += xvel*dt
    b[1] += b[4]*dt
    
    if abs(b[0]) > box_width:
      b[3] = -b[3]
      b[0] += b[3]*dt
    
    if abs(b[1]) > box_width:
      b[4] = -b[4]
      b[1] += b[4]*dt
    
    return b

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
            list2[0]+=5
    elif list1 == -1:
        if list2[0]>0:
            list2[0]-=5
    elif list1 == 10:
        if list2[1]<box_width:
            list2[1]+=5
    elif list1 == -10:
        if list2[1] > 0:
            list2[1]-=5
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
