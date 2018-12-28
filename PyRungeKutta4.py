#!/usr/bin/env python

from numpy import *
    
def RG4_D(f, X0, X1, xe, t):
    N = 10000     # Cambiar este no. de iteraciones si se requiere
    c = 2.998 * 10**(8.) # Usamos MKS
    kpc = 3.0857 * 10**19.
    kyear = 3.154 * 10**10.
    t, niter = 0, 0
    h = 10. * kpc/c/N
    D = 0
        
    Rk_a = [[1./2., 0., 0., 0.],[0., 1./2., 0., 0.], [0., 0., 1., 0.], [0., 0., 0., 0.]]
    Rk_b = [1./6., 1./3., 1./3., 1./6.]
    Rk_c = [0., 1./2., 1./2., 1./6.]
    k = zeros(4)
    k[0] = 0
    Y0 = X0

    while (D < xe * kpc):
        alpha = 0.0 * kpc
        D = sqrt((X0[0]-X1[0]) * (X0[0]-X1[0]) + (X0[1]-X1[1]) * (X0[1]-X1[1]) + (X0[2]-X1[2]) * (X0[2]-X1[2]))
        
        k1 = f( X0, t )[0]
        k2 = f( X0 + ( 0.5 * k1 * h ), t )[0]
        k3 = f( X0 + ( 0.5 * k2 * h ), t )[0]
        k4 = f( X0 + ( k3 * h ), t )[0]
        X0 = X0 + ( ( h/6.0 ) * ( k1 + ( k2 + k3 )*2. ) + k4 ) 
        t += h
        niter += 1 
        
    return X0, t, niter 