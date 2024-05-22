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

import numpy as np
import gillespie as gp
from gillespy2.solvers.cpp import SSACSolver
import scipy.spatial as spatial
from scipy.stats import kde
#import cell
import cell2 
import csv
from multiprocessing import Pool, Process, freeze_support
#from sklearn.neighbors import NearestNeighbors
import params
import time





def main_fn(core):
    np.random.seed(random.randint(1,100000))
    
	
    # mu=params.mus[core]#np.random.choice(params.mus)#[core]
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
    
    # print(myx)
 
    nd=[]
    gly_n=[]
  
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
        
        print("############ Time ################")
        
        if round(tnow,1)%1==0:
            print([tnow,len(bb),len(g3p_xy),len(gly_xy)])
            
        
        # tree=spatial.KDTree(gly_xy)
        # tree_g3p=spatial.KDTree(g3p_xy)
        # if round(tnow,1)%10==0:
        #     print(core,to_switch)
        
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
            
            
            
             
            if round(tnow,1)%5==0:
                
                            
     
    ########### uptake of glycerol from nearby   
                if params.gly_uptake==1:
                    pass
                            
                       
                         
    ########### uptake of g3p from nearby  // no osmotic
                if params.g3p_uptake==1:
                    pass
                   
               
        if tnow==0.1:
            

            fln='init_core_%s' % core
            filename1 = fln + "_histogram_%s.csv" % tnow
            
            # dat1=[[b.pmd,b.pmk,b.pmr,b.pd,b.pk,b.pr,b.pg3p,b.pg3pr,b.pgly] for b in bb]
            dat1=[[b.mu, b.pmd,b.pmr,b.pd,b.pr,b.pg3p] for b in bb]
            f1 = open(filename1, 'w')
    
            writer1 = csv.writer(f1)
    
            writer1.writerows(dat1)

            f1.close()
                    
                    
        # if round(tnow,1)%50==0:

        #     fln='tnow_core_%s' % core
        #     filename1 = fln + "_histogram_%s.csv" % tnow
            
        #     # dat1=[[b.pmd,b.pmk,b.pmr,b.pd,b.pk,b.pr,b.pg3p,b.pg3pr,b.pgly] for b in bb]
        #     dat1=[[b.pmd,b.pmr,b.pd,b.pr,b.pg3p] for b in bb]
        #     f1 = open(filename1, 'w')
    
        #     writer1 = csv.writer(f1)
    
        #     writer1.writerows(dat1)

        #     f1.close()
            
    
        
        directions = [1,-1,10,-10]
        step=random.choices(directions,k=len(g3p_xy))
            
        # g3p_xy=cell2.update(step,g3p_xy,'G3P')
        
        if len(g3p_xy)>0:
            g3p_xy=np.array(list(map(cell2.update_g3p,step,g3p_xy)))
        
        step=random.choices(directions,k=len(gly_xy))
            
        # gly_xy=cell2.update(step,gly_xy,'Gly')
        
        if len(gly_xy)>0:
            gly_xy=np.array(list(map(cell2.update_gly,step,gly_xy)))
        
    
    hist1=[[b.pd,b.pmd,b.pg3p, b.pmr, b.pr, b.ppd, b.ppr] for b in bb]
    # hist2=[b.pmd for b in bb]
    # density = kde.gaussian_kde(hist)
    # x = np.linspace(min(hist),max(hist),100)
    # y=density(x)
    
    # plt.plot(x, y)
    # plt.title("Density Plot of the data")
    # plt.show()
    # print(np.mean(hist))
    
    # gra=[b.gmax for b in bb]
    # print(gra)
  
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
    
    

    g3p_xy=np.array(g3p_xy)
    
    


if __name__ == '__main__':
    with Pool(51) as p:
        print(p.map(main_fn, range(1,51)))




