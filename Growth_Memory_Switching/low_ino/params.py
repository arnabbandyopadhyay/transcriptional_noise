#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 22, 2022

@author: Arnab Bandyopadhyay
"""

import random
import numpy as np

box_width=500
random.seed(5)

mus=[random.randrange(30,120) for rn in range(2000)]
sigma = 10
dt = 0.1
max_cells=10000
inheritance=1 # 
gmax=[random.randrange(30,40) for rn in range(200)]
myx=[random.uniform(0.1, 0.5) for i in range(1500)]
reten_prob=[random.uniform(0, 0.5) for i in range(100)]
switching_prob=[random.uniform(0, 0.5) for i in range(100)]
switching_back=[random.uniform(0, 0) for i in range(100)]
radius=3
max_lim=[5,50]
method='stochastic'
ran=np.random.uniform(0,3,1000)
factors=[1,1]#[ff for ff in ran if ff>0]
km=100

