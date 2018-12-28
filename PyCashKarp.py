#!/usr/bin/env python

from numpy import *

# Cash-Karp parameters
A = [ 0.0, 0.2, 0.3, 0.6, 1.0, 0.875 ]
B = zeros((6,6))
B[1,:] = [0.2,0,0,0,0,0]
B[2,:] = [3.0/40.0, 9.0/40.0,0,0,0,0]
B[3,:] = [0.3, -0.9, 1.2,0,0,0]
B[4,:] = [-11.0/54.0, 2.5, -70.0/27.0, 35.0/27.0,0,0]
B[5,:] = [1631.0/55296.0, 175.0/512.0, 575.0/13824.0, 44275.0/110592.0, 253.0/4096.0,0]
C  = [37.0/378.0, 0.0, 250.0/621.0, 125.0/594.0, 0.0, 512.0/1771.0]
DC = [C[0]-2825.0/27648.0, C[1]-0.0, C[2]-18575.0/48384.0,C[3]-13525.0/55296.0, C[4]-277.00/14336.0, C[5]-0.25]

def Cash_Karp(f, X0, X1, xe, t):
    c = 2.998*10**(8.) # Usamos MKS
    kpc = 3.0857*10**19.
    kyear = 3.154*10**10.
    t, niter = 0, 0 
    niter_si_avanza, niter_no_avanza = 0, 0
    
    hmax = 0.50*kpc/c # 400 pc/c
    delta_error_max = 0.0001 
    D = 0
    K = zeros((6,6))
    h = hmax
    n_trial = 0
    while (D < xe*kpc ):
        if D <= 20*kpc and n_trial == 0:
            h = hmax
        elif (D > 20*kpc and n_trial==0):
            h = hmax*1000
        
        D = sqrt((X0[0]-X1[0]) * (X0[0]-X1[0]) + (X0[1]-X1[1]) * (X0[1]-X1[1]) + (X0[2]-X1[2]) * (X0[2]-X1[2]))
        K[:,0] = h*f(X0,t)[0]
        K[:,1] = h*f(X0 + B[1,0]*K[:,0], t + h*A[1])[0]
        K[:,2] = h*f(X0 + B[2,1]*K[:,1] + B[2,0]*K[:,0], t + h*A[2])[0]
        K[:,3] = h*f(X0 + B[3,2]*K[:,2] + B[3,1]*K[:,1] + B[3,0]*K[:,0], t + h*A[3])[0]
        K[:,4] = h*f(X0 + B[4,3]*K[:,3] + B[4,2]*K[:,2] + B[4,1]*K[:,1] + B[4,0]*K[:,0], t + h*A[4])[0]
        K[:,5] = h*f(X0 + B[5,4]*K[:,4] + B[5,3]*K[:,3] + B[5,2]*K[:,2] + B[5,1]*K[:,1] + B[5,0]*K[:,0], t + h*A[5])[0]
        
        desp = (C[1]*K[0:2,0]+ C[1]*K[0:2,1]+ C[2]*K[0:2,2]+ C[3]*K[0:2,3]+ C[4]*K[0:2,4]+ C[5]*K[0:2,5]).max()
        error = (DC[0]*K[0:2,0]+ DC[1]*K[0:2,1]+ DC[2]*K[0:2,2]+ DC[3]*K[0:2,3]+ DC[4]*K[0:2,4]+ DC[5]*K[0:2,5]).max()
        delta_error = abs(error/desp)
        niter += 1
        if delta_error < delta_error_max:
            X0 = X0 + C[0]*K[:,0] + C[1]*K[:,1] + C[2]*K[:,2] + C[3]*K[:,3] + C[4]*K[:,4] + C[5]*K[:,5]
            t += h
            niter_si_avanza += 1
            h = hmax
            n_trial = 0
        else:
            h = h/2.
            niter_no_avanza += 1
            n_trial=1

    return X0, t, niter, niter_si_avanza 