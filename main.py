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
#import cell
import cell2 
import csv
from multiprocessing import Pool, Process, freeze_support

import time
# def uptake(a,b):
    
#     for i in range(len(b)):
        
#         if g3px_xy[b[i],2]>2:
#             a.pg3p += 2
#             g3px_xy[b[i],2] -= 2
#             break



random.seed(5)
mus=[random.randrange(45,120) for rn in range(100)]
sigma = 20

random.seed()
fig = plt.figure()
camera = Camera(fig)



def main_fn(core):
    
	
    mu=mus[core]
    n_particles = 100
    box_width = 500
    # n_steps = 10#500#700
    dt = 1
    tnow=0
    max_cells=10000
    radius=10
    
    # def uptake(a,b):
    
    #     for i in range(len(b)):
            
    #         if g3px_xy[b[i],2]>2:
    #             a.pg3p += 2
    #             g3px_xy[b[i],2] -= 2
    #             break
    
    nd=[]
    gly_n=[]
    # g3p_xy=[[0,0]]
    # g3px_xy=np.array([[0,0,0,0.1,0.1]])
    
    bb=[cell2.Cell(0, mu) for i in range(50)]
    
    gly_xy=np.random.randint(box_width,size=[100000,2])
    g3p_xy=np.random.randint(box_width,size=[10,2])
    # plt.plot(gly_xy[:,0],gly_xy[:,1],'r.',markersize=1)
    # plt.savefig('gly_initial.png')
    
    
    # for j in range(n_steps):
    while len(bb)<max_cells:
        
        # print(len(gly_xy))
        tnow=tnow+1
        print([tnow,len(bb),len(g3p_xy)])
        
        tree=spatial.KDTree(gly_xy)
        tree_g3p=spatial.KDTree(g3p_xy)
        
        cell2.Cell.tnow=tnow
        
        bb=[b.proliferate() if tnow > b.nextdiv else [b] for b in bb]
        
        bb=sum(bb,[])
        
        bb_xy=[[b.xpos,b.ypos] for b in bb]
        
        clpt_gly=tree.query_ball_point(bb_xy,radius)
        clpt_g3p=tree_g3p.query_ball_point(bb_xy,radius)
    
        counter=-1
        for i in bb:
            counter+=1
            i.volume=i.volume*np.exp(np.log(2)*dt/int(i.nextdiv-i.birthtime))
            # if len(clpt_gly[counter])>0:
            #     length=len(clpt_gly[counter])
            #     if length>=2:
            #         i.pgly+=2
            #         gly_xy=np.delete(gly_xy,clpt_gly[counter][:2],axis=0)
            #         tree=spatial.KDTree(gly_xy)
            #         clpt_gly=tree.query_ball_point(bb_xy,5)
            #     elif length==1:
            #         i.pgly+=1
            #         gly_xy=np.delete(gly_xy,clpt_gly[counter][:1],axis=0)
            #         tree=spatial.KDTree(gly_xy)
            #         clpt_gly=tree.query_ball_point(bb_xy,5)
            
            if tnow%2==0:
     
    ############ uptake of glycerol from nearby               
                if len(clpt_gly[counter])>0:
                    length=int(0.05*len(clpt_gly[counter])*i.pk)
                    # print('length')
                    # print(length)
                    
                    if length>=0:
                        i.pgly+=length
                        gly_xy=np.delete(gly_xy,clpt_gly[counter][:length],axis=0)
                        tree=spatial.KDTree(gly_xy)
                        clpt_gly=tree.query_ball_point(bb_xy,radius)
                    
                    
                    # if length>=5:
                    #     i.pgly+=5
                    #     gly_xy=np.delete(gly_xy,clpt_gly[counter][:5],axis=0)
                    #     tree=spatial.KDTree(gly_xy)
                    #     clpt_gly=tree.query_ball_point(bb_xy,radius)
                    # else:
                    #     i.pgly+=length
                    #     gly_xy=np.delete(gly_xy,clpt_gly[counter][:length],axis=0)
                    #     tree=spatial.KDTree(gly_xy)
                    #     clpt_gly=tree.query_ball_point(bb_xy,radius)
                        
                        
    ############## Osmotic uptake and release of g3p
                g3pn=cell2.osmotic(i.pg3p,i.volume, len(clpt_g3p[counter]), radius)
                if g3pn>i.pg3p:
                    # print('uptake')
                    in_out_diff=int(g3pn-i.pg3p)
                    g3p_xy=np.delete(g3p_xy,clpt_g3p[counter][:in_out_diff],axis=0)
                    tree_g3p=spatial.KDTree(g3p_xy)
                    clpt_g3p=tree_g3p.query_ball_point(bb_xy,radius)
                    i.pg3p=g3pn
                else:
                    in_out_diff=int(i.pg3p-g3pn)
                    # print('release')
                    # print(in_out_diff)
                    if in_out_diff>0:
                        g3p_cor=[[i.xpos,i.ypos] for b in range(in_out_diff)]
                        g3p_xy=np.append(g3p_xy,g3p_cor,axis=0)
                        tree_g3p=spatial.KDTree(g3p_xy)
                        clpt_g3p=tree_g3p.query_ball_point(bb_xy,radius)
                        i.pg3p=g3pn
                        
                
                
                kk=i.ssa()
            
                    
# ############## uptake of g3p from nearby

#             if len(clpt_g3p[counter])>0:
#                 length=len(clpt_g3p[counter])
#                 # print('length')
#                 # print(length)
#                 if length>=5:
#                     i.pg3p+=5
#                     g3p_xy=np.delete(g3p_xy,clpt_g3p[counter][:5],axis=0)
#                     tree_g3p=spatial.KDTree(g3p_xy)
#                     clpt_g3p=tree_g3p.query_ball_point(bb_xy,radius)
#                 else:
#                     i.pg3p+=length
#                     g3p_xy=np.delete(g3p_xy,clpt_g3p[counter][:length],axis=0)
#                     tree_g3p=spatial.KDTree(g3p_xy)
#                     clpt_g3p=tree_g3p.query_ball_point(bb_xy,radius)
                    
                    
#             g3pn=round(0.1*np.random.random()*i.pg3p)
#             # print(g3pn)
#             if g3pn>0:
#                 g3p_cor=[[i.xpos,i.ypos] for b in range(g3pn)]
#                 g3p_xy=np.append(g3p_xy,g3p_cor,axis=0)
#                 i.pg3p-=g3pn
            
#             # print(g3p_xy)
            
            # if tnow%1==0:
            #     kk=i.ssa()
            
            
                # g3px_xy=np.append(g3px_xy, [[i.xpos,i.ypos, kk,
                #                              2*(np.random.random()-0.5),
                #                              2*(np.random.random()-0.5)]],axis=0)
                
        
        
            
            
            
            
            # print([i.xpos,i.ypos])
            # print([[x,y] for x,y in zip(g3px_x,g3px_y)])
            
################ Velocity based update
            
            # i.xpos += i.xvel*dt
            # i.ypos += i.yvel*dt
        
            # if abs(i.xpos) > box_width:
            #   i.xvel = -i.xvel
            #   i.xpos += i.xvel*dt
              
            # if i.xpos < 0:
            #   i.xvel = -i.xvel
            #   i.xpos += i.xvel*dt
            
            
            # if abs(i.ypos) > box_width:
            #   i.yvel = -i.yvel
            #   i.ypos += i.yvel*dt
              
            # if i.ypos < 0:
            #   i.yvel = -i.yvel
            #   i.ypos += i.yvel*dt
              
################# Random walk update cells

            directions = ["UP", "DOWN", "LEFT", "RIGHT"]
            step = random.choice(directions)
            
            if step == "RIGHT":
                if i.xpos<box_width:
                    i.xpos+=5
            elif step == "LEFT":
                if i.xpos>0:
                    i.xpos-=5
            elif step == "UP":
                if i.ypos<box_width:
                    i.ypos+=5
            elif step == "DOWN":
                if i.ypos > 0:
                    i.ypos-=5
                    
                    
################# Random walk update g3p

        # start = time.time()
        
        
        directions = [1,-1,10,-10]
        step=random.choices(directions,k=len(g3p_xy))
            
        # g3p_xy=cell2.update(step,g3p_xy,'G3P')
        
        g3p_xy=list(map(cell2.update_g3p,step,g3p_xy))
        
        step=random.choices(directions,k=len(gly_xy))
            
        # gly_xy=cell2.update(step,gly_xy,'Gly')
        
        gly_xy=list(map(cell2.update_gly,step,gly_xy))
        
        # print(time.time()-start)



        
        # if j==0:
        #     g3px_xy=np.delete(g3px_xy,0,axis=0)
        
        # g3px_xy=np.array(g3px_xy)
        
        # g3px_xy=np.array(list(map(cell.take_step,g3px_xy)))
        
        # point_tree = spatial.cKDTree(g3px_xy[:,[0,1]])
        
        # clpt=[point_tree.query_ball_point([b.xpos,b.ypos], 1) for b in bb]
        
        # list(map(uptake,bb,clpt))
        
        x_coord=[b.xpos for b in bb]
        y_coord=[b.ypos for b in bb]
        
        nd.append([x_coord,y_coord])
        
        # print([b.pg3p for b in bb])
        # print(len(g3p_xy))
        
        # if j%10 == 0:
        #     plt.plot(nd[j][0],nd[j][1],'ro')
        #     camera.snap()
            
    #     if j%10 == 0:
    #         xcord=[b[0] for b in g3px_xy]
    #         ycord=[b[1] for b in g3px_xy]
    #         plt.plot(nd[j][0],nd[j][1],'ro')
    #         plt.plot(xcord, ycord,'b.')
    #         camera.snap()
    
    #     if j%10 == 0:
    #         xcord=[b[0] for b in g3p_xy]
    #         ycord=[b[1] for b in g3p_xy]
    #         plt.plot(xcord, ycord,'b.')
    #         camera.snap()
            
    # animation = camera.animate()
    # animation.save('celluloid_minimal_2.gif', writer = 'imagemagick') 
          
           
            # print([i.xpos,i.ypos,i.xvel,i.yvel])
    
    # Creating histogram
    
    hist1=[[b.pd,b.pmd,b.pg3p] for b in bb]
    # hist2=[b.pmd for b in bb]
    # density = kde.gaussian_kde(hist)
    # x = np.linspace(min(hist),max(hist),100)
    # y=density(x)
    
    # plt.plot(x, y)
    # plt.title("Density Plot of the data")
    # plt.show()
    # print(np.mean(hist))
    
    filename1 = "histogram_pd_%s.csv" % core
    # filename2 = "histogram_pmd_%s.csv" % core
    # filename = "histogram_1.csv"
    f1 = open(filename1, 'w')
    # f2 = open(filename2, 'w')
    
    # create the csv writer
    writer1 = csv.writer(f1)
    # writer2 = csv.writer(f2)
    
    # write a row to the csv file
    writer1.writerows(hist1)
    # writer2.writerow(hist2)
    
    # close the file
    f1.close()
    # f2.close()
    
    
    file = "gly_%s.csv" % core
    f = open(file, 'w')
    wr = csv.writer(f)
    nm=[len(gly_xy)]
    wr.writerow(nm)
    f.close()
    
    
    
    
    
    plt.plot(gly_xy[:,0],gly_xy[:,1],'r.',markersize=1)
    
    file = "gly_%s.png" % core
    plt.savefig(file)
    plt.close()
    
    
    
    plt.plot(g3p_xy[:,0],g3p_xy[:,1],'r.',markersize=1)
    file = "g3p_%s.png" % core
    plt.savefig(file)
    plt.close()
    
    print(len(g3p_xy))
   
    

# main_fn(1)

if __name__ == '__main__':
    with Pool(10) as p:
        print(p.map(main_fn, range(1,11)))



# if __name__ == '__main__':
#     freeze_support()
#     p1 = Process(target=main_fn, args=(1,))
#     p1.start()
#     p1.join()

#     p2 = Process(target=main_fn, args=(2,))
#     p2.start()
#     p2.join()
