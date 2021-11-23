# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 16:29:15 2021

@author: Arnab
"""
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib import animation
import numpy as np
import random
from celluloid import Camera

mu, sigma = 7, 3 

class Cell:
    
    def __init__(self, birthtime, **kwargs):
        
        
        self.birthtime=birthtime
        
         
        self.nextdiv=self.birthtime+np.random.normal(mu, sigma)#np.random.normal(mu, sigma)
        
        # self.generation = 0
        
        if 'xpos' in kwargs: self.xpos = kwargs['xpos']
        else: self.xpos = round(np.random.random(),2)*box_width
        
        if 'ypos' in kwargs: self.ypos = kwargs['ypos']
        else: self.ypos =round(np.random.random(),2)*box_width
        
        if 'xvel' in kwargs: self.xvel = kwargs['xvel']
        else: self.xvel = 2*(np.random.random()-0.5)*box_width
        
        if 'yvel' in kwargs: self.yvel = kwargs['yvel']
        else: self.yvel =2*(np.random.random()-0.5)*box_width
        
            
        
    @classmethod  
    def clone(cls, b):
        return cls(birthtime=tnow, nextdiv=tnow+np.random.normal(mu, sigma), xpos = b.xpos,
                   ypos = b.ypos,xvel=2*(np.random.random()-0.5)*box_width,yvel=2*(np.random.random()-0.5)*box_width)
    def gen(self):
        self.generation += 1

    def proliferate(self):


        nl=[Cell.clone(self),Cell.clone(self)]

        return nl


# def get_initial_coordinates():
#   x_coord = [np.random.random()*box_width for i in range(n_particles)]
#   y_coord = [np.random.random()*box_width for i in range(n_particles)]
  
#   return x_coord, y_coord


# def get_initial_velocities():
    
#     # x_vel=[250 for i in range(n_particles)]
#     # y_vel=[250 for i in range(n_particles)]
#     x_vel = [2*(np.random.random()-0.5)*box_width for i in range(n_particles)]
#     y_vel = [2*(np.random.random()-0.5)*box_width for i in range(n_particles)]
#     return x_vel, y_vel

def take_step(x_coord, y_coord, x_vel, y_vel):
  for i in range(n_particles):
    x_coord[i] += x_vel[i]*dt
    y_coord[i] += y_vel[i]*dt
    
    if abs(x_coord[i]) > box_width:
      x_vel[i] = -x_vel[i]
      x_coord[i] += x_vel[i]*dt
    
    if abs(y_coord[i]) > box_width:
      y_vel[i] = -y_vel[i]
      y_coord[i] += y_vel[i]*dt
    
  return x_coord, y_coord, x_vel, y_vel

# def take_step(b):
    
#     b.xpos += b.xpos*dt
#     b.ypos += b.ypos*dt
        
#     if abs(b.xpos) > box_width:
#         b.xvel = -b.xvel
#         b.xpos += b.xvel*dt
      
#     if abs(b.ypos) > box_width:
#         b.yvel = -b.yvel
#         b.ypos += b.yvel*dt
  
#     return b        
            
    
    
    
  
      
      
      

# def add_frame(xs,ys,i):
#   global trajectory
#   if i == 0:
#     trajectory = ''
#   trajectory += str(n_particles) + '\ntitle\n'
#   for x, y in zip(xs,ys):
#     trajectory += ' '.join(['Ar',str(x),str(y),'0.0\n'])
    
    
import numpy as np

n_particles = 100
box_width = 10
n_steps = 250
dt = 0.1
tnow=0
global trajectory

# x_coord, y_coord = get_initial_coordinates()

# x_vel, y_vel = get_initial_velocities()

fig = plt.figure()
camera = Camera(fig)

bb=[Cell(0),Cell(0)]

nd=[]
nd2=[]
for j in range(n_steps):
    print(j)
    
    tnow=j*dt
    
    ll=[b.proliferate() if tnow > b.nextdiv else [b] for b in bb]
    
    ll=sum(ll,[])
    
    # n_particles=len(bb)
    
    for i in ll:
        i.xpos += i.xvel*dt
        i.ypos += i.yvel*dt
    
        if abs(i.xpos) > box_width:
          i.xvel = -i.xvel
          i.xpos += i.xvel*dt
        else:
          i.xvel = i.xvel
        
        if abs(i.ypos) > box_width:
          i.yvel = -i.yvel
          i.ypos += i.yvel*dt
        else:
          i.yvel = i.yvel
    
    x_coord=[b.xpos for b in ll]
    y_coord=[b.ypos for b in ll]
    x_vel=[b.xvel for b in ll]
    # y_vel=[b.yvel for b in ll]
    print(x_vel)
    # x_coord, y_coord, x_vel, y_vel = take_step(x_coord, y_coord, x_vel, y_vel)
    # print(x_coord)
    # ll= [take_step(b) for b in ll]
    
    # zz=[b.xpos for b in ll]
    # zz2=[b.ypos for b in ll]
    # nd.append([zz,zz2,tnow])
    
    # zz=list(x_coord)
    # zz2=list(y_coord)
    nd.append(x_coord)
    nd2.append(y_coord)
    
    # if j%10 == 0:
    #     add_frame(x_coord, y_coord,j)
    
    if j%10 == 0:
        plt.plot(nd[j],nd2[j],'ro')
        camera.snap()
        
animation = camera.animate()
animation.save('celluloid_minimal.gif', writer = 'imagemagick') 
      
      
     
        
    
# import py3Dmol
# view = py3Dmol.view()
# view.addModelsAsFrames(trajectory,'xyz')
# view.animate({'loop': 'forward', 'reps': 1})
# view.setStyle({'sphere':{'radius': 0.1}})
# view.zoomTo()

# n=4
# particles = np.zeros(n, dtype=[("position", float, 2)])

# particles["position"] = [[i,j]for i,j in zip(nd[0][0],nd[0][1])]

# particles["size"] = 0.5 * np.ones(n);


# fig = plt.figure(figsize=(10, 10))
# ax = plt.axes(xlim=(-10, 12), ylim=(-10, 12))
# scatter = ax.scatter(particles["position"][:, 0], particles["position"][:, 1],s=100,c='#ff7f0e')



# def update(ii):
#     particles["position"] = [[i,j]for i,j in zip(nd[ii][0],nd[ii][1])]
   
#     scatter.set_offsets(particles["position"])
#     return scatter,

# Writer = animation.writers['ffmpeg']
# writer = Writer(fps=5, metadata=dict(artist='Me'), bitrate=1800)
# anim = FuncAnimation(fig, update, interval=10)
# anim.save('gc.mp4',writer=writer)
# plt.show()


# ll=[Cell(0)]
# for i in range(1000):
    
#     ll=[take_step(b) for b in ll]
#     print([[b.xpos, b.ypos, b.xvel, b.yvel] for b in ll])

