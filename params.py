#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 22 09:03:46 2022

@author: abp19
"""

import random
import numpy as np


box_width=500
random.seed(5)

mus=[30,90,30,90,30,90,30,90,30,90,30,90,30,90]#[random.randrange(30,91,60) for rn in range(20)]#[random.randrange(30,90) for rn in range(100)]
mus=[random.randrange(30,120) for rn in range(2000)]
sigma = 10
dt = 0.1
max_cells=1000
inheritance=1 # 
gmax=[random.randrange(30,35) for rn in range(200)]
myx=[random.uniform(0.001, 0.5) for i in range(500)]#(0.001, 0.5)(0.3, 0.5)
radius=3
gly_uptake=0 #from environment=1, 0 otherwise
g3p_uptake=0 #from environment=1, 0 otherwise
max_lim=[5,50]
method='stochastic'
# ran=np.random.normal(3,3,1000)
ran=np.random.uniform(0,3,1000)
factors=[1,1]#[ff for ff in ran if ff>0]
