# -*- coding: utf-8 -*-
"""
Created on Thu Nov 25 23:42:05 2021

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
import cell 
import csv
from multiprocessing import Pool, Process, freeze_support


# def uptake(a,b):
    
#     for i in range(len(b)):
        
#         if g3px_xy[b[i],2]>2:
#             a.pg3p += 2
#             g3px_xy[b[i],2] -= 2
#             break




mu, sigma = 120, 3

# fig = plt.figure()
# camera = Camera(fig)



def main_fn(core):
    

    n_particles = 100
    box_width = 1000
    n_steps = 900
    dt = 1
    tnow=0
    
    
    
    def uptake(a,b):
    
        for i in range(len(b)):
            
            if g3px_xy[b[i],2]>2:
                a.pg3p += 2
                g3px_xy[b[i],2] -= 2
                break
    
    nd=[]
    
    g3px_xy=np.array([[0,0,0,0.1,0.1]])
    
    bb=[cell.Cell(0) for i in range(5)]
    
    for j in range(n_steps):
        print(j)
        tnow=dt*j
        cell.Cell.tnow=tnow
        
        bb=[b.proliferate() if tnow > b.nextdiv else [b] for b in bb]
        
        bb=sum(bb,[])
        
        
    
    
        
        for i in bb:
            kk=i.ssa()
            
            
            g3px_xy=np.append(g3px_xy, [[i.xpos,i.ypos, kk,
                                         2*(np.random.random()-0.5),
                                         2*(np.random.random()-0.5)]],axis=0)
            
            
            # print([i.xpos,i.ypos])
            # print([[x,y] for x,y in zip(g3px_x,g3px_y)])
            
            i.xpos += i.xvel*dt
            i.ypos += i.yvel*dt
        
            if abs(i.xpos) > box_width:
              i.xvel = -i.xvel
              i.xpos += i.xvel*dt
            
            
            if abs(i.ypos) > box_width:
              i.yvel = -i.yvel
              i.ypos += i.yvel*dt
        
        if j==0:
            g3px_xy=np.delete(g3px_xy,0,axis=0)
        
        # g3px_xy=np.array(g3px_xy)
        
        g3px_xy=np.array(list(map(cell.take_step,g3px_xy)))
        
        point_tree = spatial.cKDTree(g3px_xy[:,[0,1]])
        
        clpt=[point_tree.query_ball_point([b.xpos,b.ypos], 1) for b in bb]
        
        list(map(uptake,bb,clpt))
        
        x_coord=[b.xpos for b in bb]
        y_coord=[b.ypos for b in bb]
        
        nd.append([x_coord,y_coord])
        
        # if j%10 == 0:
        #     plt.plot(nd[j][0],nd[j][1],'ro')
        #     camera.snap()
            
    #     if j%10 == 0:
    #         xcord=[b[0] for b in g3px_xy]
    #         ycord=[b[1] for b in g3px_xy]
    #         plt.plot(nd[j][0],nd[j][1],'ro')
    #         plt.plot(xcord, ycord,'b.')
    #         camera.snap()
            
    # animation = camera.animate()
    # animation.save('celluloid_minimal_2.gif', writer = 'imagemagick') 
          
           
            # print([i.xpos,i.ypos,i.xvel,i.yvel])
    
    # Creating histogram
    
    hist=[b.pd for b in bb]
    # density = kde.gaussian_kde(hist)
    # x = np.linspace(min(hist),max(hist),100)
    # y=density(x)
    
    # plt.plot(x, y)
    # plt.title("Density Plot of the data")
    # plt.show()
    # print(np.mean(hist))
    
    filename = "histogram_%s.csv" % core
    # filename = "histogram_1.csv"
    f = open(filename, 'w')
    
    # create the csv writer
    writer = csv.writer(f)
    
    # write a row to the csv file
    writer.writerow(hist)
    
    # close the file
    f.close()




if __name__ == '__main__':
    with Pool(2) as p:
        print(p.map(main_fn, [1, 2, 3]))



# if __name__ == '__main__':
#     freeze_support()
#     p1 = Process(target=main_fn, args=(1,))
#     p1.start()
#     p1.join()

#     p2 = Process(target=main_fn, args=(2,))
#     p2.start()
#     p2.join()
