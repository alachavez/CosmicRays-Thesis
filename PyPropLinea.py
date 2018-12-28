#!/usr/bin/env python

from numpy import *
from PyPropagador import *

data = loadtxt("../Datos/Datos_StarBurst.txt");
  
Energyi = data[2,2]; loni = data[2,5]; lati = -data[2,6];

galactic_l_rad = deg2rad(loni);
galactic_b_rad = deg2rad(lati);
nx = 30*cos( galactic_b_rad ) * cos( galactic_l_rad )
ny = 30*cos( galactic_b_rad ) * sin( galactic_l_rad ) 
nz = 30*sin( galactic_b_rad )
    

#Propagador(Charge, Mass, Energy, galactic_l, galactic_b, xf, yf, zf, BField, Distance)
Energyf, lonf, latf, xf, yf, zf = Propagador(1, 1, 80, loni, lati, nx, ny, nz, "BField", 30)

print(Energyf, lonf, latf, xf, yf, zf)
