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
from sklearn.neighbors import NearestNeighbors
import params
import time

# def uptake(a,b):
    
#     for i in range(len(b)):
        
#         if g3px_xy[b[i],2]>2:
#             a.pg3p += 2
#             g3px_xy[b[i],2] -= 2
#             break





# random.seed()
# fig = plt.figure()
# camera = Camera(fig)



def main_fn(core):
    np.random.seed(random.randint(1,100000))
    
	
    mu=params.mus[core]#np.random.choice(params.mus)#[core]
    myx=np.random.choice(params.myx)
    sigma=params.sigma
    method=params.method
    
    box_width = params.box_width
    # n_steps = 10#500#700
    dt = params.dt
    tnow=0
    max_cells=params.max_cells
    radius=params.radius
    lmt=np.random.choice(params.max_lim)
    to_switch=np.random.binomial(1,1)
    prob=np.random.choice(params.switching_prob)
    sw_back=np.random.choice(params.switching_back)
    reten_prob=np.random.choice(params.reten_prob)
    print(myx)
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
#    if np.random.uniform(0,1)>0.5:
#        bb=[cell2.Cell(0, np.random.choice(params.gmax), lmt, np.random.choice(params.factors), np.random.choice(params.myx1)) for i in range(100)]
#    else:
#        bb=[cell2.Cell(0, np.random.choice(params.gmax), lmt, np.random.choice(params.factors), np.random.choice(params.myx2)) for i in range(100)]
#    
    bb=[cell2.Cell(0, np.random.choice(params.gmax), lmt, np.random.choice(params.factors), np.random.choice(params.myx)) for i in range(20)]
    
    if params.gly_uptake==1: 
        gly_xy=np.random.randint(box_width,size=[200000,2])
    else:
        gly_xy=np.random.randint(box_width,size=[5,2])
    
    
    if params.g3p_uptake==1: 
        g3p_xy=np.random.randint(box_width,size=[75000,2])
    else:
        g3p_xy=np.random.randint(box_width,size=[5,2])
    # plt.plot(gly_xy[:,0],gly_xy[:,1],'r.',markersize=1)
    # plt.savefig('gly_initial.png')
    
    
    # for j in range(n_steps):
    while len(bb)<max_cells:
        
        # print(len(gly_xy))
        tnow=round(tnow+dt,1)
        
        
        
        if round(tnow,1)%1==0:
            print([tnow,len(bb),len(g3p_xy),len(gly_xy)])
            
        
        # tree=spatial.KDTree(gly_xy)
        # tree_g3p=spatial.KDTree(g3p_xy)
        if round(tnow,1)%10==0:
            print(core,to_switch)
        
        cell2.Cell.tnow=tnow
        
        if ((to_switch==1) and (len(bb)<max_cells)):
            switch='yes'
        else:
            switch='no'
        
        bb=[b.proliferate(switch, prob, sw_back, reten_prob) if tnow > b.nextdiv else [b] for b in bb]
        
        
        
        
        bb=sum(bb,[])
        
#        bbl=[b.influx for b in bb if b.influx <=0.3]
#        bbh=[b.influx for b in bb if b.influx >0.3]
#        print([len(bbl),len(bbh)])
        
        bb_xy=[[b.xpos,b.ypos] for b in bb]
        
        
        
        # if len(gly_xy)>0:
        #     neigh_gly = NearestNeighbors(radius=radius)
        #     neigh_gly.fit(gly_xy)
        #     clpt_gly=np.asarray(neigh_gly.radius_neighbors(bb_xy,return_distance=True,sort_results=True)[1])
            
        # if len(g3p_xy)>0:
        #     neigh_g3p = NearestNeighbors(radius=radius)
        #     neigh_g3p.fit(g3p_xy)
        #     clpt_g3p=np.asarray(neigh_g3p.radius_neighbors(bb_xy,return_distance=True,sort_results=True)[1])
        
        
        # clpt_gly=neigh_gly.kneighbors(bb_xy,return_distance=False)
        # clpt_gly=np.asarray(neigh_gly.radius_neighbors(bb_xy,return_distance=True,sort_results=True)[1])
        # clpt_g3p=np.asarray(neigh_g3p.radius_neighbors(bb_xy,return_distance=True,sort_results=True)[1])
        # clpt_g3p=neigh_g3p.kneighbors(bb_xy, return_distance=False)
        # clpt_gly=tree.query_ball_point(bb_xy,radius)
        # clpt_g3p=tree_g3p.query_ball_point(bb_xy,radius)
        
        
        
        # if round(tnow,1)%5==0:
        #     if method=='deterministic':
        #            [i.deterministic_ode() for i in bb]
        
        
    
        counter=-1
        for i in bb:
#            if i.mu<30:
#                print(['new born mu ',i.mu])
            
            counter+=1
            i.volume=1
            #i.volume=i.volume+dt*(np.log(2)/int(i.nextdiv-i.birthtime))*np.exp(np.log(2)*(tnow-i.birthtime-1)/int(i.nextdiv-i.birthtime))
            
            ################# Random walk update cells

            directions = ["UP", "DOWN", "LEFT", "RIGHT"]
            step = random.choice(directions)
            
            if step == "RIGHT":
                if i.xpos<box_width:
                    i.xpos+=2
            elif step == "LEFT":
                if i.xpos>0:
                    i.xpos-=2
            elif step == "UP":
                if i.ypos<box_width:
                    i.ypos+=2
            elif step == "DOWN":
                if i.ypos > 0:
                    i.ypos-=2
            
            
            
            
            
            
            # if round(tnow,1)%1==0:
            #     print([i.volume,i.nextdiv])
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
            
            if round(tnow,1)%5==0:
                
    # ######################### osmotic uptake glycerol
    #             glyn=cell2.osmotic(i.pgly,i.volume, len(clpt_gly[counter]), radius)
    #             # print([i.pgly,len(clpt_gly[counter]),glyn])
    #             if glyn>i.pgly:
    #                 # print('uptake')
    #                 in_out_diff=int(glyn-i.pgly)
    #                 gly_xy=np.delete(gly_xy,clpt_gly[counter][:in_out_diff],axis=0)
    #                 # neigh_gly = NearestNeighbors(radius=radius)
    #                 # neigh_gly.fit(gly_xy)
    #                 # clpt_gly=neigh_gly.kneighbors(bb_xy, return_distance=False)
    #                 tree=spatial.KDTree(gly_xy)
    #                 clpt_gly=tree.query_ball_point(bb_xy,radius)
    #                 # print('uptake')
    #                 # print([i.pgly,len(clpt_gly[counter]),glyn])
    #                 # print(in_out_diff)
    #                 i.pgly=glyn
                            
     
    ########### uptake of glycerol from nearby   
                if params.gly_uptake==1:
                    pass
                            
                #     if len(clpt_gly[counter])>0:
                #         if i.pk>0:
                            
                #             length=int(1*len(clpt_gly[counter])*np.random.normal(i.pk,0.1))
                            
                #             if length>len(clpt_gly[counter]):
                #                 length=len(clpt_gly[counter])
                                
                                
                #             # print('length')
                #             # print(length)
                            
                #             if length>=0:
                #                 i.pgly+=length
                #                 gly_xy=np.delete(gly_xy,clpt_gly[counter][:length],axis=0)
                                
                #                 if len(gly_xy)>0:
                #                     neigh_gly = NearestNeighbors(radius=radius)
                #                     neigh_gly.fit(gly_xy)
                #                     clpt_gly=np.asarray(neigh_gly.radius_neighbors(bb_xy,return_distance=True,sort_results=True)[1])
                                    
                                
                                
    
                                # clpt_gly=neigh_gly.kneighbors(bb_xy, return_distance=False)
                                # tree=spatial.KDTree(gly_xy)
                                # clpt_gly=tree.query_ball_point(bb_xy,radius)
                            
                         
    ########### uptake of g3p from nearby  // no osmotic
                if params.g3p_uptake==1:
                    pass
                    # length_g3p=len(clpt_g3p[counter])
                    # i.pg3p+=length_g3p
                    # g3p_xy=np.delete(g3p_xy,clpt_g3p[counter][:length_g3p],axis=0)
                    # if len(g3p_xy)>0:
                    #     neigh_g3p = NearestNeighbors(radius=radius)
                    #     neigh_g3p.fit(g3p_xy)
                    #     clpt_g3p=np.asarray(neigh_g3p.radius_neighbors(bb_xy,return_distance=True,sort_results=True)[1])
                    
                    # range_x=g3p_xy[(g3p_xy[:,0]>=i.xpos-radius) & (g3p_xy[:,0]<=i.xpos+radius)]
                    # range_y=range_x[(range_x[:,1]>=i.ypos-radius) & (range_x[:,1]<=i.ypos+radius)]
                    
                    # to_add=int(np.random.uniform(0,0.3,1)*len(np.unique(range_y,axis=0))) #int((i.pt/(5+i.pt))*0.5*len(np.unique(range_y,axis=0)))
                    # print(to_add)
                    
                    # i.pg3p+=to_add
                    # t_delete=range_y[np.random.choice(range(len(np.unique(range_y,axis=0))),to_add)]
                    # print(range_y)
                    # g3p_xy=np.delete(g3p_xy,range_y[0:to_add],axis=0)
                
                    # if len(clpt_g3p[counter])>0:
                    #     if i.pk>0:
                    #         length_g3p=int(1*len(clpt_g3p[counter])*np.random.normal(i.pk,0.1))
                    #         if length_g3p>len(clpt_g3p[counter]):
                    #             length_g3p=len(clpt_g3p[counter])
                                
                    #         if length_g3p>=0:
                    #             i.pg3p+=length_g3p
                    #             g3p_xy=np.delete(g3p_xy,clpt_g3p[counter][:length_g3p],axis=0)
                                
                    #             if len(g3p_xy)>0:
                    #                 neigh_g3p = NearestNeighbors(radius=radius)
                    #                 neigh_g3p.fit(g3p_xy)
                    #                 clpt_g3p=np.asarray(neigh_g3p.radius_neighbors(bb_xy,return_distance=True,sort_results=True)[1])
        
                                
                                
                    
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
                        
                        
    # ############## Osmotic uptake and release of g3p
    #             g3pn=cell2.osmotic(i.pg3p,i.volume, len(clpt_g3p[counter]), radius)
    #             if g3pn>i.pg3p:
    #                 # print('uptake')
    #                 in_out_diff=int(g3pn-i.pg3p)
    #                 g3p_xy=np.delete(g3p_xy,clpt_g3p[counter][:in_out_diff],axis=0)
    #                 neigh_g3p = NearestNeighbors(radius=radius)
    #                 neigh_g3p.fit(g3p_xy)
    #                 clpt_g3p=np.asarray(neigh_g3p.radius_neighbors(bb_xy,return_distance=True,sort_results=True)[1])
    #                 # clpt_g3p=neigh_g3p.kneighbors(bb_xy, return_distance=False)
    #                 # tree_g3p=spatial.KDTree(g3p_xy)
    #                 # clpt_g3p=tree_g3p.query_ball_point(bb_xy,radius)
    #                 i.pg3p=g3pn
                # else:
                #     in_out_diff=int(i.pg3p-g3pn)
                #     # print('release')
                #     # print(in_out_diff)
                #     if in_out_diff>0:
                #         g3p_cor=[[i.xpos,i.ypos] for b in range(in_out_diff)]
                #         g3p_xy=np.append(g3p_xy,g3p_cor,axis=0)
                #         neigh_g3p = NearestNeighbors(radius=radius)
                #         neigh_g3p.fit(g3p_xy)
                #         clpt_g3p=np.asarray(neigh_g3p.radius_neighbors(bb_xy,return_distance=True,sort_results=True)[1])
                #         # clpt_g3p=neigh_g3p.kneighbors(bb_xy, return_distance=False)
                #         # tree_g3p=spatial.KDTree(g3p_xy)
                #         # clpt_g3p=tree_g3p.query_ball_point(bb_xy,radius)
                #         i.pg3p=g3pn
                        
                
                # if method=='stochastic':
                #     kk=i.ssa()
                #     # print([i.pmd,i.pmr,i.pd,i.pr,i.pg3p])
                # if method=='deterministic':
                #     kk=i.simple_ode()
                # if method=='sde':
                #     kk=i.simple_sde()
                    
                    

                    
                
        # if len(bb)>max_cells-500:
    
        #     fln='core_%s' % core
        #     filename1 = fln + "_histogram_pd_%s.csv" % tnow
            
        #     # dat1=[[b.pmd,b.pmk,b.pmr,b.pd,b.pk,b.pr,b.pg3p,b.pg3pr,b.pgly] for b in bb]
        #     dat1=[[b.pmd,b.pmr,b.pd,b.pr,b.pg3p] for b in bb]
        #     f1 = open(filename1, 'w')
    
        #     writer1 = csv.writer(f1)
    
        #     writer1.writerows(dat1)

        #     f1.close()
                    
               
        if tnow==0.1:
            

            fln='init_core_%s' % core
            filename1 = fln + "_histogram_%s.csv" % tnow
            
            # dat1=[[b.pmd,b.pmk,b.pmr,b.pd,b.pk,b.pr,b.pg3p,b.pg3pr,b.pgly] for b in bb]
            dat1=[[b.mu, b.pmd,b.pmr,b.pd,b.pr,b.pg3p] for b in bb]
            f1 = open(filename1, 'w')
    
            writer1 = csv.writer(f1)
    
            writer1.writerows(dat1)

            f1.close()
                    
                    
        if round(tnow,1)%50==0:

            fln='tnow_core_%s' % core
            filename1 = fln + "_histogram_%s.csv" % tnow
            
            # dat1=[[b.pmd,b.pmk,b.pmr,b.pd,b.pk,b.pr,b.pg3p,b.pg3pr,b.pgly] for b in bb]
            dat1=[[b.pmd,b.pmr,b.pd,b.pr,b.pg3p] for b in bb]
            f1 = open(filename1, 'w')
    
            writer1 = csv.writer(f1)
    
            writer1.writerows(dat1)

            f1.close()
            
            
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
              

                    
################# Random walk update g3p

        # start = time.time()
        
        
        directions = [1,-1,10,-10]
        step=random.choices(directions,k=len(g3p_xy))
            
        # g3p_xy=cell2.update(step,g3p_xy,'G3P')
        
        if len(g3p_xy)>0:
            g3p_xy=np.array(list(map(cell2.update_g3p,step,g3p_xy)))
        
        step=random.choices(directions,k=len(gly_xy))
            
        # gly_xy=cell2.update(step,gly_xy,'Gly')
        
        if len(gly_xy)>0:
            gly_xy=np.array(list(map(cell2.update_gly,step,gly_xy)))
        
        # print(time.time()-start)
        
        

        
        # if j==0:
        #     g3px_xy=np.delete(g3px_xy,0,axis=0)
        
        # g3px_xy=np.array(g3px_xy)
        
        # g3px_xy=np.array(list(map(cell.take_step,g3px_xy)))
        
        # point_tree = spatial.cKDTree(g3px_xy[:,[0,1]])
        
        # clpt=[point_tree.query_ball_point([b.xpos,b.ypos], 1) for b in bb]
        
        # list(map(uptake,bb,clpt))
        
        # x_coord=[b.xpos for b in bb]
        # y_coord=[b.ypos for b in bb]
        
        # nd.append([x_coord,y_coord])
        
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
    
    hist1=[[b.pd,b.pmd,b.pg3p, b.pmr, b.pr, b.ppd, b.ppr] for b in bb]
    # hist2=[b.pmd for b in bb]
    # density = kde.gaussian_kde(hist)
    # x = np.linspace(min(hist),max(hist),100)
    # y=density(x)
    
    # plt.plot(x, y)
    # plt.title("Density Plot of the data")
    # plt.show()
    # print(np.mean(hist))
    
  
    filename1 = "histogram_pd_%s.csv" % core
    print('entered')
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
    
    
    # file = "gly_%s.csv" % core
    # f = open(file, 'w')
    # wr = csv.writer(f)
    # nm=[len(gly_xy)]
    # wr.writerow(nm)
    # f.close()
    
    g3p_xy=np.array(g3p_xy)
    
    
    # plt.plot(gly_xy[:,0],gly_xy[:,1],'r.',markersize=1)
    
    # file = "gly_%s.png" % core
    # plt.savefig(file)
    # plt.close()
    
    
    
#    plt.plot(g3p_xy[:,0],g3p_xy[:,1],'r.',markersize=1)
#    file = "g3p_%s.png" % core
#    plt.savefig(file)
#    plt.close()
    
    
   
    

# main_fn(1)

if __name__ == '__main__':
    with Pool(10) as p:
        print(p.map(main_fn, range(31,51)))



# if __name__ == '__main__':
#     freeze_support()
#     p1 = Process(target=main_fn, args=(1,))
#     p1.start()
#     p1.join()

#     p2 = Process(target=main_fn, args=(2,))
#     p2.start()
#     p2.join()

