#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 22 09:03:46 2022

@author: abp19
"""

import random

box_width=500
random.seed(5)
mus=[random.randrange(30,35) for rn in range(100)]
sigma = 5
dt = 0.1
max_cells=10000
radius=3
gly_uptake=1 #from environment=1, 0 otherwise
g3p_uptake=1 #from environment=1, 0 otherwise
max_lim=[5,25]
method='deterministic'