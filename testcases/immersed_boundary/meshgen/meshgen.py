#! /usr/bin/env python

import numpy as np
import scipy as sp
import os.path
import iris.cube
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import code
from meshfuncs import *



fnm_in = 'ib_config.mcf'

fin = open(fnm_in,'r')
lines = fin.readlines()
fin.close()

for l in lines:
  if 'x_size' in l:nx = np.int(l.strip('\n').split("=")[-1])
  if 'y_size' in l:ny = np.int(l.strip('\n').split("=")[-1])
  if 'dxx' in l:dx = np.float(l.strip('\n').split("=")[-1])
  if 'dyy' in l:dy = np.float(l.strip('\n').split("=")[-1])




#-----------------------------------------
# Random 2D periodic field
if False:
#if True:
  lamb = 300.      # wavelength in m
  hmax = 500.    # peak height in m
  hmin = 5.    # minumim height in m
  f = 2
  halo = 2
  seed = -1 # set to -1 for no seed
  nodes, elems, norms = Random2DOperator(f*nx,f*ny,dx/f,dy/f, lamb, hmax, hmin, halo, seed)
  

#-----------------------------------------
# cosine hill (Wood 95)
#if False:
if True:
  lamb = 500.      # wavelength in m
  hmax = 500.    # peak height in m
  dim_type = 2   # 1=1D, 2=2D
  f = 2
  offsetx = 0.
  offsety = 0.
  offsetz = 0.
  cx = dx*np.float(nx)/2.0
  cy = dy*np.float(ny)/2.0
  noise = False
  halo = 2
  nodes, elems, norms = Wood95Operator(f*nx,f*ny,dx/f,dy/f, cx, cy, lamb, hmax, dim_type, offsetx, offsety, offsetz, noise, halo)
  

#-----------------------------------------
# blocks
if False:
#if True:
  cents = [] # center location
  hdims = []  # x & y dimensions
  hghts = []  # heights
  dxs = dx/2.
  dys = dy/2.
  halo = 2
  offx = dxs/2.
  offy = dys/2.
  offz = 5. # min height
  cents.append((dx*np.float(nx)/3. , dy*np.float(ny)/3.))
  cents.append((2*dx*np.float(nx)/3. , 2*dy*np.float(ny)/3.))
  hdims.append((150., 150.)) # this will be snapped to nearest grid node
  hdims.append((200., 200.)) # this will be snapped to nearest grid node
  hghts.append(500.)
  hghts.append(200.)
  nodes, elems, norms = BlocksOperator(cents, hdims, hghts, dxs, dys, dx*np.float(nx), dy*np.float(ny), offx,offy, offz, halo)
   

#-----------------------------------------
## Flat surface for testing
if False:
#if True:
  hmax = 5.2
  f=2
  offsetx = 0
  offsety = 0
  halo = 2
  nodes, elems, norms = FlatOperator(hmax, f*nx,f*ny,dx/f,dy/f, offsetx, offsety, halo)
  

#-----------------------------------------
# Bolund data
#if True:
if False:
  dpath = '/net/home/h02/tdunstan/monc/meshgen/Bolund/Bolund_data.npz' 
  dxs = 0.25    # source data dx (m)
  dys = 0.25    # source data dy (m)
  f = 2 # sample spacing frequency for output data (dx_out=dx/f)
  rot=31 # rotation angle for datai
  halo = 2 # must be >= MONC halo size
  xoff = 0.0 # x offset
  yoff = 0.0 # y offset
  zoff = 0.7 # z offset (minimum height)
  gsig = (dx/f)/dxs # gaussian filter sigma in grid lengths (0.0 = off)
  print('sigma =',gsig)
  nodes, elems, norms = BolundRotOperator(dpath, f*nx, f*ny, dx/f, dy/f, halo, xoff, yoff, zoff, dxs, dys, rot, gsig)



#########################################################################################

cubes = []
nnodes = nodes.shape[1]
nelems = elems.shape[1]
cubes.append(iris.cube.Cube(nnodes, var_name='nnodes'))
cubes.append(iris.cube.Cube(nelems, var_name='nelems'))
cubes.append(iris.cube.Cube(nodes, var_name='nodes'))
cubes.append(iris.cube.Cube(elems, var_name='elems'))
cubes.append(iris.cube.Cube(norms, var_name='norms'))
iris.save(cubes,'mesh.nc')

# plotting
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
tmpx = [] 
tmpy = [] 
tmpz = [] 
for i in range(nnodes):
  if nodes[2,i] > -0.1:
    tmpx.append(nodes[0,i])
    tmpy.append(nodes[1,i])
    tmpz.append(nodes[2,i])
ax.scatter(tmpx, tmpy, tmpz, c='r', s=10., edgecolors='none')
plt.show()
#

