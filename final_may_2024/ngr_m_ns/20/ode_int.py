#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  8 09:03:49 2022

@author: abp19
"""

import numpy as np
from scipy.integrate import odeint,solve_ivp
import matplotlib.pyplot as plt



def simple_model(t,z):
    
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
    
    
    
    
    kon_pD_R = 2165.44314226246
    kon_pT_R = 1588.06318455853
    koff_pD_R = 10014.7094615103 
    koff_pT_R = 10477.3841801672
    
    ktrsc = 0.5 
    r_R = 2
    ktrsl = 0.3 
    ktrsp_G3P = 0 # 0.01757360326999
        
    kon_G3P_R = 99.7805018542808
                        
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
    return dzdt



# # # initial condition
# z0 = [0.5,0.5,100,100,100,100,100,100,100]

# # # # # number of time points
# n = 1000
# tf=1000
# t = np.linspace(0,tf,n)




# z = solve_ivp(simple_model,(0,tf),z0,t_eval=np.linspace(0,tf,n).tolist(),method='LSODA')
    

# # plot results
# plt.plot(t,z.y[5,:],'g:',label='u(t)')
# plt.ylabel('values')
# plt.xlabel('time')
# plt.legend(loc='best')
# plt.show()







# function that returns dz/dt
# def model(t,z):
    
#     pD_R=z[0]
#     pFK_R=z[1] 
#     pT_R=z[2]
#     mD=z[3]
#     mF=z[4]
#     mK=z[5]
#     mT=z[6] 
#     mR=z[7] 
#     D=z[8]
#     F=z[9]
#     K=z[10]
#     T=z[11]
#     R=z[12]
#     Gly=z[13]
#     G3P=z[14]
#     G3PR=z[15]
#     GlyK=z[16]
#     G3PD=z[17]
    
#     # pD_R, pFK_R, pT_R, mD, mF, mK, mT, mR, D, F, K, T, R, Gly, G3P, G3PR, GlyK, G3PD=z
    
    
    
    
#     kon_pD_R = 12.0999
#     kon_pFK_R = 11.0653  
#     kon_pT_R = 7.10524
#     koff_pD_R = 86.7934 
#     koff_pFK_R = 79.1819 
#     koff_pT_R = 37.8405
#     ktrsp_G3P = 0 #0.0100652           
#     k_glycoly = 0.0153946           
#     kon_G3P_R = 11.4709              
#     kon_Gly_K = 37.1818              
#     kon_G3P_D = 0.01               
#     koff_G3P_R = 0.360667            
#     koff_Gly_K = 0.794284           
#     koff_G3P_D = 0.278724            
#     kc_GlyK = 7.50307               
#     kc_G3PD = 82.4072
#     ktrsc = 0.5                           
#     r_a_i = 20                                           
#     r_R = 2                                              
#     ktrsl = 0.3 
#     dRNA = 0.069                 
#     dPROTEIN = 0.023                
#     dG3P =  0.035                                        
#     dG3PR = 0.035                
#     pR = 1                                
#     pD_t= 1                                 
#     pFK_t = 1                             
#     pT_t = 1                                                                                   
#     pD = pD_t - pD_R                                     
#     pFK= pFK_t- pFK_R                                    
#     pT = pT_t - pT_R    
            
#     ktrsp_Gly = 0 # 0.01 
    
    
    
    
    
#     ra1 = kon_pD_R * R *(pD_t - pD_R) - koff_pD_R * pD_R
#     ra2 = kon_pFK_R* R *(pFK_t- pFK_R)- koff_pFK_R* pFK_R
#     ra3 = kon_pT_R * R *(pT_t - pT_R) - koff_pT_R * pT_R
#     rb1 = ktrsc * pD 
#     rb2 = ktrsc * pFK
#     rb3 = ktrsc * pT 
#     rb4 = ktrsc / r_a_i * pD_R 
#     rb5 = ktrsc / r_a_i * pFK_R
#     rb6 = ktrsc / r_a_i * pT_R
#     rb7 = ktrsc / r_R * pR
#     rc1 = ktrsl * mD 
#     rc2 = ktrsl * mF
#     rc3 = ktrsl * mK
#     rc4 = ktrsl * mT 
#     rc5 = ktrsl * mR
#     rd1 = ktrsp_Gly * F
#     rd2 = ktrsp_G3P * T + k_glycoly
#     re1 = kon_G3P_R * G3P * R - koff_G3P_R * G3PR
#     re2 = kon_Gly_K * Gly * K - koff_Gly_K * GlyK
#     re3 = kc_GlyK * GlyK
#     re4 = kon_G3P_D * G3P * D - koff_G3P_D * G3PD
#     re5 = kc_G3PD * G3PD 
#     rf1 = dRNA * mD
#     rf2 = dRNA * mF
#     rf3 = dRNA * mK
#     rf4 = dRNA * mT
#     rf5 = dRNA * mR
#     rf6 = dPROTEIN * D
#     rf7 = dPROTEIN * F
#     rf8 = dPROTEIN * K
#     rf9 = dPROTEIN * T 
#     rf10= dPROTEIN * R 
#     rf11= dG3P * G3P
#     rf12= dG3PR* G3PR 
                     

    
#     pD_Rt  = ra1
#     pFK_Rt = ra2
#     pT_Rt  = ra3
#     mDt = rb1 + rb4 - rf1
#     mFt = rb2 + rb5 - rf2  
#     mKt = rb2 + rb5 - rf3  
#     mTt = rb3 + rb6 - rf4 
#     mRt = rb7 - rf5 
#     Dt = rc1 - re4 + re5 - rf6
#     Ft = rc2 - rf7 
#     Kt = rc3 - re2 + re3 - rf8    
#     Tt = rc4 - rf9
#     Rt = rc5 - re1 - rf10
#     Glyt  = rd1 - re2 
#     G3Pt  = rd2 - re1 + re3 - re4 - rf11
#     G3PRt = re1 - rf12
#     GlyKt = re2 - re3 
#     G3PDt = re4 - re5               

                            

#     dzdt = [pD_Rt, pFK_Rt, pT_Rt, mDt, mFt, mKt, mTt, mRt, Dt, Ft, Kt,
#             Tt, Rt, Glyt, G3Pt, G3PRt, GlyKt, G3PDt]
#     return dzdt

# # # initial condition
# z0 = [0,0,0,0,5,5,5,5,5,5,5,5,5,5,5,5,5,5]

# # # # # number of time points
# # n = 401

# # t = np.linspace(0,40,n)




# z = solve_ivp(model,(0,5),z0,t_eval=np.linspace(0,5,100).tolist())
    

# # plot results
# plt.plot(t,z[:,18],'g:',label='u(t)')
# plt.plot(t,x,'b-',label='x(t)')
# plt.plot(t,y,'r--',label='y(t)')
# plt.ylabel('values')
# plt.xlabel('time')
# plt.legend(loc='best')
# plt.show()
