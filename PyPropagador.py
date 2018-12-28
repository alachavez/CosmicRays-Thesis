#!/usr/bin/env python

from numpy import *
from PyRungeKutta4 import RG4_D
from PyCashKarp import Cash_Karp
from PyJansonFarrar import *


def Propagador(Charge, Mass, Energy, galactic_l, galactic_b, xf, yf, zf, BField, Distance):
    
    
    #if BField == "JanssonFarrar": getField = JanssonFarrar.getField
    #else if: BField == "Pshirkov"; getField = Pshirkov.getField
    #else if: BField == "Sun2008"; getField = Sun2008.getField  
    #else if: BField == "Stanev1"; getField = Stanev1.getField
    #else if: BField == "Stanev2"; getField = Stanev2.getField
    #else if: BField == "Stanev3"; getField = Stanev3.getField


    # Entradas:  Charge,Mass,Energy, galactic_l, galactic_b in degrees

    EeV = 1.602 * 10**(-1.)
    Mp = 1.67 * 10**(-27.)
    Qp = 1.602 * 10**(-19.) 
    c = 2.998 * 10**(8.)
    kpc = 3.0857 * 10**19.
    Mpc = kpc * 10**3.
    muGauss = 10.**(-10.) # Tesla
    nGauss = muGauss * 10**(-3.)
    kyear = 3.154 * 10**10.
    
    #if CM == 'sup'
      #  kpc = Mpc;
        #  muGauss = nGauss;
    #end
 
    #----- Cond. iniciales -----
    E0 = Energy * EeV  # p=Ev/c^2
    q = Charge * Qp   # units of p charge, negative foe backpropagation
    m0 = Mass * Mp     # units of p mass
     
    galactic_l_rad = deg2rad( galactic_l )
    galactic_b_rad = deg2rad( galactic_b )
    
    nx = cos( galactic_b_rad ) * cos( galactic_l_rad )
    ny = cos( galactic_b_rad ) * sin( galactic_l_rad ) 
    nz = sin( galactic_b_rad )
    n_inicial = [nx, ny, nz]
    n_inicial = n_inicial/linalg.norm( n_inicial )
    
    v0 = c * sqrt( 1.0 - m0**2 * c**4 / E0**2) * n_inicial
    p_i = E0 * v0 / c**2
    x0 = xf * kpc #-8.5 * kpc;
    y0 = yf * kpc
    z0 = zf * kpc
    px0 = p_i[0] 
    py0 = p_i[1] 
    pz0 = p_i[2]
    X0_Rel = [x0, y0, z0, px0, py0, pz0]
 
    def f_Rel(X, t):
        M = m0 * sqrt( 1 + (X[3]**2 + X[4]**2 + X[5]**2) / (m0 * c)**2); 
        x = X[0] / kpc
        y = X[1] / kpc
        z = X[2] / kpc
        (Bx, By, Bz) = getField(x, y, z) 
        px = X[3] / M
        py = X[4] / M
        pz = X[5] / M
        ax = q/M*(X[4]*Bz* muGauss - X[5]*By* muGauss)
        ay = q/M*(X[5]*Bx* muGauss - X[3]*Bz* muGauss)
        az = q/M*(X[3]*By* muGauss - X[4]*Bx* muGauss)
        return [px, py, pz, ax, ay, az]
    
    t0 = 0.
    Xe = Distance
    Xi_Rel = [x0, y0, z0]
    xt, t_f, niter, niter2 = Cash_Karp(f_Rel, X0_Rel, Xi_Rel, Xe, t0)
    x_f = xt[0] / kpc 
    y_f = xt[1] / kpc 
    z_f = xt[2] / kpc
    px = xt[3] 
    py = xt[4] 
    pz = xt[5]
    p_f = [px, py, pz]
    
    M = m0 * sqrt( 1 + (px**2 + py**2 + pz**2) / (m0 * c)**2)
    E_f = M * c**2 / EeV
    t_f = t_f #/ kyear;
    R_f = sqrt(x_f**2 + y_f**2 + z_f**2)

    if dot(n_inicial, p_f) / linalg.norm( p_f ) == 1.0: 
        alpha = rad2deg(arccos(dot(n_inicial, p_f) / linalg.norm( p_f ) ) )
    else:
        alpha = rad2deg(arccos(0.99) )
    
    l_f = rad2deg( arctan2( py, px ) )
    if l_f < 0: l_f += 360.
    b_f = rad2deg( arctan2( pz, sqrt(px * px + py * py ) ) )
    
    return ([E_f, l_f, b_f, x_f, y_f, z_f])
    #return ([E_f])
