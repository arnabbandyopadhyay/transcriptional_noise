# -*- coding: utf-8 -*-
"""
Created on Thu Nov 25, 2021

@author: Arnab Bandyopadhyay
"""

import numpy as np
import random
import cell
import csv
from multiprocessing import Pool, Process, freeze_support
import params
import time


def main_fn(core):
    np.random.seed(random.randint(1,100000))
    
    sigma=params.sigma
    method=params.method
    
    box_width = params.box_width
    dt = params.dt
    tnow=0
    max_cells=params.max_cells
    lmt=np.random.choice(params.max_lim)
    to_switch=np.random.binomial(1,1)
    prob=np.random.choice(params.switching_prob)
    sw_back=np.random.choice(params.switching_back)
    reten_prob=np.random.choice(params.reten_prob)

 
    nd=[]
    gly_n=[]
  
    bb=[cell.Cell(0, np.random.choice(params.gmax), lmt, np.random.choice(params.factors), np.random.choice(params.myx)) for i in range(100)]
    
    gly_xy=np.random.randint(box_width,size=[5,2])
    g3p_xy=np.random.randint(box_width,size=[5,2])
    
    while len(bb)<max_cells:
        
        tnow=round(tnow+dt,1)
        
        print("############ Time ################")
        
        if round(tnow,1)%1==0:
            print([tnow,len(bb),len(g3p_xy),len(gly_xy)])
        
        cell.Cell.tnow=tnow
        
        if ((to_switch==1) and (len(bb)<max_cells)):
            switch='yes'
        else:
            switch='no'
        
        bb=[b.proliferate(switch, prob, sw_back, reten_prob) if tnow > b.nextdiv else [b] for b in bb]
        
        bb=sum(bb,[])
        
        bb_xy=[[b.xpos,b.ypos] for b in bb]
    
        counter=-1
        for i in bb:
            
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
                   
               
        if tnow==0.1:
            
            fln='init_core_%s' % core
            filename1 = fln + "_histogram_%s.csv" % tnow
            
            dat1=[[b.mu, b.pmd,b.pmr,b.pd,b.pr,b.pg3p] for b in bb]
            f1 = open(filename1, 'w')
    
            writer1 = csv.writer(f1)
    
            writer1.writerows(dat1)

            f1.close()
        
    
    hist1=[[b.pd,b.pmd,b.pg3p, b.pmr, b.pr, b.ppd, b.ppr] for b in bb]
  
    filename1 = "histogram_pd_%s.csv" % core
    print('computation finished, writing in files')
    f1 = open(filename1, 'w')
    writer1 = csv.writer(f1)
    writer1.writerows(hist1)
    f1.close()
    
    


if __name__ == '__main__':
    with Pool(21) as p:
        print(p.map(main_fn, range(1,21)))




