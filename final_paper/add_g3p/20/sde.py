#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 15:19:25 2022
https://pypi.org/project/sdeint/
@author: abp19
"""

# Zombie apocalypse SDE model

import matplotlib.pyplot as plt
import numpy as np
import sdeint

# P, d, B, G, A = 0.0001, 0.0001, 0.0095, 0.0001, 0.0001

# tspan = np.linspace(0, 2., 1000)
# y0 = np.array([10., 10., 10., P])


# def f(y, t):
#     Si = y[0]
#     Zi = y[1]
#     Ri = y[2]

#     f0 = y[3] - B * Si * Zi - d * Si
#     f1 = B * Si * Zi + G * Ri - A * Si * Zi
#     f2 = d * Si + A * Si * Zi - G * Ri
#     f3 = 0
#     return np.array([f0, f1, f2, f3])


# def GG(y, t):
#     return np.diag([1, 1, 1, 1])

# result = sdeint.itoint(f, GG, y0, tspan)


# plt.plot(result)
# plt.show()

def simple_model(z,t):
    
    pD_R=z[0]
    pT_R=z[1]
    mD=z[2]
    mT=z[3] 
    mR=z[4] 
    D=z[5]
    T=z[6]
    R=z[7]
    G3P=z[8]
    
    # pD_R, pFK_R, pT_R, mD, mF, mK, mT, mR, D, F, K, T, R, Gly, G3P, G3PR, GlyK, G3PD=z
    
    
    scaling=1e9
    
    kon_pD_R = 2165.44314226246/scaling
    kon_pT_R = 1588.06318455853/scaling
    koff_pD_R = 10014.7094615103/scaling
    koff_pT_R = 10477.3841801672/scaling
    
    ktrsc = 0.5 
    r_R = 2
    ktrsl = 0.3 
    ktrsp_G3P = 0.01757360326999
        
    kon_G3P_R = 99.7805018542808/scaling
                        
    kon_G3P_D = 0.103423441900033

    dRNA = 0.069              
    dPROTEIN = 0.023                
    dG3P =  0.035    
    pR = 1                                
    pD_t= 1                                 
    pT_t = 1                                                                                   
    pD = pD_t - pD_R                                                                  
    pT = pT_t - pT_R    
            
    
    ra1 = kon_pD_R * R *(pD_t - pD_R) - koff_pD_R * pD_R
    ra3 = kon_pT_R * R *(pT_t - pT_R) - koff_pT_R * pT_R    
    rb1 = ktrsc * pD     
    rb3 = ktrsc * pT 
    rb7 = ktrsc / r_R * pR
    rc1 = ktrsl * mD 
    rc4 = ktrsl * mT 
    rc5 = ktrsl * mR
    rd2 = ktrsp_G3P * T
    re1 = kon_G3P_R * G3P * R
    re4 = kon_G3P_D * G3P * D 
    rf1 = dRNA * mD
    rf4 = dRNA * mT
    rf5 = dRNA * mR
    rf6 = dPROTEIN * D
    rf9 = dPROTEIN * T 
    rf10= dPROTEIN * R 
    rf11= dG3P * G3P
    
                     

    
    pD_Rt  = ra1
    pT_Rt  = ra3 
    mDt = rb1 - rf1
    mTt = rb3 - rf4 
    mRt = rb7 - rf5 
    Dt = rc1 - rf6    
    Tt = rc4 - rf9
    Rt = rc5 - re1 - rf10    
    G3Pt  = rd2 - re1 - re4 - rf11             

                            

    dzdt = [pD_Rt, pT_Rt, mDt, mTt, mRt, Dt,
            Tt, Rt, G3Pt]
    return np.array(dzdt)


def GG(y, t):
    return np.diag([0, 0, 1, 1, 1, 1, 1, 1, 1])

# result = sdeint.itoEuler(simple_model, GG, [0.5,0.5,10.,10.,10.,10.,10.,10.,10.], tspan)


# plt.plot(result)
# plt.show()