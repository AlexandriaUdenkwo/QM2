#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 19:46:36 2023

@author: alexandriaudenkwo
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as spy
import scipy.special as special

h = 0.001
E = 1
m = 939 #MeV
l = 1
hbar = spy.value('reduced Planck constant times c in MeV fm') #MeV


Vo = -51 #MeV
ro = 1.27 #fm
a = 0.67 #fm
A = 16 #for oxygen
Ro = ro*(A**(1/3))

#intial conditions
u_0 = 0
u_1 = 10**(-14)


def V(r):
     return Vo/(1 + np.exp((r - Ro)/a)) 

def g(r):
    if 0 in r:
        for i in range(len(r)):
            if r[i]==0:
                r[i] = h/2
        return ((2*m)/(hbar**2))*(E - V(r) - (l*(l+1)*(hbar**2))/(2*m*((r)**2)))
    else:
        return ((2*m)/(hbar**2))*(E - V(r) - (l*(l+1)*(hbar**2))/(2*m*((r)**2)))
 
#phase shifts points
r1, r2 = 10, 15
       
    
def main():
    r = np.linspace(-5, 5, 10001)
    u = np.zeros(r.size)
    u[0], u[1] = u_0, u_1
    g_r = g(r)
    #print(g_r)

    for i in range(len(r)-2):   
   
        u[i+2] = (2*u[i+1]*(1 - (5/12)*(h**2)*(g_r[i+1])) - u[i]*(1 + (1/12)*(h**2)*(g_r[i])))/(1 + (1/12)*(h**2)*(g_r[i+2]))
    
    # plt.figure(1)
    # plt.plot(r, u)
    
    # plt.figure(2)
    # plt.plot(r, V(r))
    # plt.show()
    
    #phase shifts
    r_phase = np.linspace(-15, 15, 30001)
    u_phase = np.zeros(r_phase.size)
    u_phase[0], u_phase[1] = u_0, u_1
    u_l = np.zeros(r.size)
    u_r = np.zeros(r.size)
    g_rp = g(r_phase)
    
    gamma = np.zeros(r.size)
    
    r_l= np.zeros(r.size)
    r_r = np.zeros(r.size)
    
    k_l = np.zeros(r.size)
    k_r = np.zeros(r.size)
    
    for i in range(len(r_phase) -2):
        u_phase[i+2] = (2*u_phase[i+1]*(1 - (5/12)*(h**2)*(g_rp[i+1])) - u_phase[i]*(1 + (1/12)*(h**2)*(g_rp[i])))/(1 + (1/12)*(h**2)*(g_rp[i+2]))
    # for i in range(len(r_phase)):
    #     u_outside_r[0] = u_phase[len(r_phase)-len(r)]
    #     r_outside_l[-1] = r_phase[len(r) -1]
    #     r_outside_r[0] = r_phase[len(r_phase)-len(r)]
    #     if r_phase[i] not in r:
            
    #         # if i <= len(r)-1:
    #         #     r_outside_l[i] = r_phase[i]
    #         #     u_outside_l[i] = u_phase[i]
    #         #     #if i == len(r) -1:
    #         #         #print(i)
            
    #         if i >= len(r_phase) - len(r) -1:
    #             u_outside_r[i-len(r_phase) +len(r)] = u_phase[i]
    #             r_outside_r[i-len(r_phase) +len(r)] = r_phase[i]
    
    
    for i in range(len(r)):
        r_l[i] = r_phase[i+len(r_phase) - len(r)]
        u_l[i] = u_phase[i+len(r_phase) - len(r)]
        k_l[i] = g_rp[i+len(r_phase) - len(r)]

    
    for j in range(len(r)-1, -1, -1):
        r_r[-j-1] = r_phase[j+len(r_phase) - len(r)]
        u_r[-j-1] = u_phase[j+len(r_phase) - len(r)]
        k_r[-j-1] = g_rp[j+len(r_phase) - len(r)]
    
    #print(r_l, r_r)

            
    #print(u_outside_l, u_outside_r,max(u_outside_l), min(u_outside_r))            
                #if i == len(r_phase) - len(r) +1 :
                    #print(i)
    #print(len(u_outside_l), len(u_outside_r))
    
    #print(u_outside_l[:3],u_outside_l[-3:])
    #print(u_phase[19998:20002])
 
    
    #put in units of femtometers
    # u_outside_l = u_outside_l
    # u_outside_r = u_outside_r
    
   
    
    
    for i in range(len(r)):
        # if u_outside_r[i] == 0:
        #     print(i)
        # if u_outside_l[i] == 0:
        #     print(i)
        gamma[i] = (r_r[i]*u_l[i])/(u_r[i]*r_l[i])  
    # max_ul = max(u_outside_l)
    # #print(max(gamma))
    # index_ul = np.where(u_outside_l == max_ul)[0][0]
    # ur_nm = u_outside_r[index_ul]
    
    # rl_in = r_outside_l[index_ul]
    # rr_in = r_outside_r[index_ul]
    
    delta = np.zeros(len(r)//2)
   # print(special.spherical_jn(l,100))
    for i in range(len(gamma)//2):
        delta[i] = (special.spherical_jn(l, k_l[i]*r_l[i]) - gamma[i]*special.spherical_jn(l, k_r[i]*r_r[i]))/\
            (special.spherical_yn(l, k_l[i]*r_l[i]) - gamma[i]*special.spherical_yn(l, k_r[i]*r_r[i]))

    angle = np.arctan(delta)
    print(max(angle))
    newr = r_l[:len(r_l)//2]
    newk = k_l[:len(k_l)//2]
    kr = np.multiply(newr, newk)
    
    r_ph = np.linspace(-15,15,10001)
    scaler =r_ph*(10**(3))
    checkr = ((l)*(l+1))/(scaler**2)
    plt.figure(3)
    plt.plot(newk,angle)
    
    v_r = (2*m/(hbar**2))*V(r_ph)
    plt.figure(4)
    plt.plot(r_ph,v_r)
    
    plt.figure(5)
    plt.plot(r_ph, checkr)
    plt.show()
    
    print(max(checkr))
    print(min(v_r))

    
main()