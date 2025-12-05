#! /usr/bin/env python

import numpy as np
import scipy as sp
from scipy import signal as sg
import os.path
import iris
import code
from scipy.interpolate import griddata
import matplotlib.pyplot as plt


def isclose(a, b, rel_tol=1e-04, abs_tol=0.0):
  return abs(a[:]-b) <= max(rel_tol * max(abs(max(a[:])), abs(b)), abs_tol)


def BolundRotOperator(dpath, nxin, nyin, dx, dy, halo, xoff_in, yoff_in, zoff, dxs, dys, ang, gsig):
  from scipy.ndimage import rotate
  from scipy.ndimage import gaussian_filter

  nx = nxin+(2*halo)
  ny = nyin+(2*halo)
  field = np.zeros((nx,ny))
  dhdx = np.zeros((nx,ny))
  dhdy = np.zeros((nx,ny))
  nodes = np.zeros((3,(nx*ny)))
  elems = np.zeros((4,(nx*ny)), dtype=int)
  norms = np.zeros((3,(nx*ny)))
  in_arr = np.zeros((nx,ny), dtype=int)

  cen_orig=[98.*4., 132.*4.] # centre point or orig in grid point units
  npzfile = np.load(dpath)
  orog = npzfile['orog']
  npzfile.close()

#  plt.subplot(1,2,1)
#  plt.contourf(orog,50)
#  plt.colorbar()

  zmin0 = np.min(orog)
  zmax0 = np.max(orog)
  zmax0=10*zmax0
  dh0=zmax0-zmin0

  # apply Gaussian filter
  if gsig > 1.e-6:
    orog = gaussian_filter(orog, gsig)
    # rescale to original min and max
    zmin1 = np.min(orog)
    zmax1 = np.max(orog)
    dh1=zmax1-zmin1
    orog[:,:] = zmin0 + (orog[:,:] - zmin1)*(dh0/dh1)


#  plt.subplot(1,2,2)
#  plt.contourf(orog,50)
#  plt.colorbar()
#  plt.show()

  # extend grid as needed for rotation
  # x-dir
  d1=float(orog.shape[1])-cen_orig[0]
  fillv = orog[0,0]
  if d1>cen_orig[0]:
    ngpx = int(2.*d1)-orog.shape[1] # +ve = add to XL
    orog = np.insert(orog, 0, np.full((ngpx,orog.shape[0]),fillv), axis=1)
  else:
    ngpx = -(int(2.*d1)-orog.shape[1]) # -ve = add to XR
    orog = np.insert(orog, -1, np.full((ngpx,orog.shape[0]),fillv), axis=1)
  
  # y-dir
  d2=float(orog.shape[0])-cen_orig[1]
  if d2>cen_orig[1]:
    ngpy = int(2.*d2)-orog.shape[0] # +ve = add to YL
    orog = np.insert(orog, 0, np.full((ngpy,orog.shape[1]),fillv), axis=0)
  else:
    ngpy = -(int(2.*d2)-orog.shape[0]) # -ve = add to XR
    orog = np.insert(orog, -1, np.full((ngpy,orog.shape[1]),fillv), axis=0)
 
  # add vertical offset
  orog[:,:] = orog[:,:]+zoff
  fillv = orog[0,0]
  orog = rotate(orog, ang, reshape=False, mode='constant', cval=fillv)
 

  # target grid
  xstart = 0.0
  xstop = dx*(nxin-1)
  xcount = complex(0,nxin)
  ystart = 0.0
  ystop = dy*(nyin-1)
  ycount = complex(0,nyin)
  tgrid_x, tgrid_y = np.mgrid[xstart:xstop:xcount, ystart:ystop:ycount] # start, stop, count

  # source grid
  xoff = yoff_in+ (((xstop-xstart) - (dxs*(orog.shape[0]-1)))/2.) # centre on target grid
  yoff = xoff_in+ (((ystop-ystart) - (dys*(orog.shape[1]-1)))/2.)
  xstart = xoff + 0.0
  xstop = xstart + dxs*(orog.shape[0]-1)
  xcount = complex(0,(orog.shape[0]))
  ystart = yoff + 0.0
  ystop = ystart + dys*(orog.shape[1]-1)
  ycount = complex(0,(orog.shape[1]))
  sgrid_x, sgrid_y = np.mgrid[xstart:xstop:xcount, ystart:ystop:ycount] # start, stop, count
  points = np.zeros((sgrid_x.size,2))
  points[:,1] = sgrid_x.flatten()
  points[:,0] = sgrid_y.flatten()
  

  tmp2 = griddata(points, orog.flatten(), (tgrid_x, tgrid_y), fill_value=orog[0,0], method='linear')

  # fill halos with wrapped data
  # middle
  field[halo:-halo,halo:-halo] = tmp2[:,:]
  # edges
  field[halo:-halo,:halo] = tmp2[:,-halo:]
  field[halo:-halo,-halo:] = tmp2[:,:halo]
  field[:halo,halo:-halo] = tmp2[-halo:,:]
  field[-halo:,halo:-halo] = tmp2[:halo,:]
  # corners
  field[:halo,:halo] = tmp2[-halo:,-halo:]
  field[:halo,-halo:] = tmp2[-halo:,:halo]
  field[-halo:,-halo:] = tmp2[:halo,:halo]
  field[-halo:,:halo] = tmp2[:halo,-halo:]

  inode = -1
  for i in range(nx):
    for j in range(ny):
      inode += 1
      in_arr[i,j]=inode+1 # Fortran indexing

      # gradients
  for j in range(ny):
    dhdx[:-1,j] = np.diff(field[:,j])/dx
    dhdx[-1,j] = (field[0,j]-field[-1,j])/dx
  for i in range(nx):
    dhdy[i,:-1] = np.diff(field[i,:])/dy
    dhdy[i,-1] = (field[i,0]-field[i,-1])/dy

  inode = -1
  ielem = -1
  mag = np.sqrt((dhdx[:,:]*dhdx[:,:]) + (dhdy[:,:]*dhdy[:,:]) + 1.0)
  for i in range(nx):
    for j in range(ny):
      inode += 1
      ielem += 1
      ip1 = i+1
      jp1 = j+1
      if ip1 > nx-1: ip1 = i
      if jp1 > ny-1: jp1 = j
      nodes[0,inode]=(i-halo)*dx
      nodes[1,inode]=(j-halo)*dy
      nodes[2,inode]=field[i,j]
      elems[0,ielem]=in_arr[i,j]     # associated nodes (clockwise)
      elems[1,ielem]=in_arr[ip1,j]
      elems[2,ielem]=in_arr[ip1,jp1]
      elems[3,ielem]=in_arr[i,jp1]
#      mag = np.sqrt((dhdx[i,j]*dhdx[i,j]) + (dhdy[i,j]*dhdy[i,j]) + 1.0)
      norms[0,ielem]=-dhdx[i,j]/mag[i,j] # unit normal i
      norms[1,ielem]=-dhdy[i,j]/mag[i,j] # unit normal j
      norms[2,ielem]=1.0/mag[i,j]        # unit normal k

#  plt.subplot(2,3,1)
#  plt.contourf(field,50)
#  plt.colorbar()
#  plt.subplot(2,3,2)
#  plt.contourf(dhdx,50)
#  plt.colorbar()
#  plt.subplot(2,3,3)
#  plt.contourf(dhdy,50)
#  plt.colorbar()
#  plt.subplot(2,3,4)
#  plt.contourf(-dhdx[:,:]/mag[:,:],50)
#  plt.colorbar()
#  plt.subplot(2,3,5)
#  plt.contourf(-dhdy[:,:]/mag[:,:],50)
#  plt.colorbar()
#  plt.subplot(2,3,6)
#  plt.contourf(1.0/mag[:,:],50)
#  plt.colorbar()
#  plt.show()

  return nodes, elems, norms


def Random2DOperator(nxin, nyin, dx, dy, lamb, hmax, hmin, halo, seed):
  nx = nxin+(2*halo)
  ny = nyin+(2*halo)
  field = np.zeros((nx,ny))
  dhdx = np.zeros((nx,ny))
  dhdy = np.zeros((nx,ny))
  nodes = np.zeros((3,(nx*ny)))
  elems = np.zeros((4,(nx*ny)), dtype=int)
  norms = np.zeros((3,(nx*ny)))
  in_arr = np.zeros((nx,ny), dtype=int)

  sizex = lamb/dx
  sizey = lamb/dy
  x, y = np.mgrid[-sizex:sizex+1, -sizey:sizey+1]
  g = np.exp(-0.333*(x**2/float(sizex)+y**2/float(sizey)))
  filter = g/g.sum()

  if seed >=0:np.random.seed(seed)
  tmp = np.random.rand(nxin,nyin)
  tmp2 = sg.convolve2d(tmp, filter, mode='same', boundary='wrap', fillvalue=0)  

  tmp = (tmp2 - tmp2.min())
  tmp2 = (tmp*(hmax-hmin)/tmp.max()) + hmin

  # fill halos with wrapped data
  # middle
  field[halo:-halo,halo:-halo] = tmp2[:,:]
  # edges
  field[halo:-halo,:halo] = tmp2[:,-halo:]
  field[halo:-halo,-halo:] = tmp2[:,:halo]
  field[:halo,halo:-halo] = tmp2[-halo:,:]
  field[-halo:,halo:-halo] = tmp2[:halo,:]
  # corners
  field[:halo,:halo] = tmp2[-halo:,-halo:]
  field[:halo,-halo:] = tmp2[-halo:,:halo]
  field[-halo:,-halo:] = tmp2[:halo,:halo]
  field[-halo:,:halo] = tmp2[:halo,-halo:]

  inode = -1
  for i in range(nx):
    for j in range(ny):
      inode += 1
      in_arr[i,j]=inode+1 # Fortran indexing

      # gradients
  for j in range(ny):
    dhdx[:-1,j] = np.diff(field[:,j])/dx
    dhdx[-1,j] = (field[0,j]-field[-1,j])/dx
  for i in range(nx):
    dhdy[i,:-1] = np.diff(field[i,:])/dy
    dhdy[i,-1] = (field[i,0]-field[i,-1])/dy

  inode = -1
  ielem = -1
  for i in range(nx):
    for j in range(ny):
      inode += 1
      ielem += 1
      ip1 = i+1
      jp1 = j+1
      if ip1 > nx-1: ip1 = i
      if jp1 > ny-1: jp1 = j
      nodes[0,inode]=(i-halo)*dx
      nodes[1,inode]=(j-halo)*dy
      nodes[2,inode]=field[i,j]
      elems[0,ielem]=in_arr[i,j]     # associated nodes (clockwise)
      elems[1,ielem]=in_arr[ip1,j]
      elems[2,ielem]=in_arr[ip1,jp1]
      elems[3,ielem]=in_arr[i,jp1]
      mag = np.sqrt((dhdx[i,j]*dhdx[i,j]) + (dhdy[i,j]*dhdy[i,j]) + 1.0)
      norms[0,ielem]=-dhdx[i,j]/mag # unit normal i
      norms[1,ielem]=-dhdy[i,j]/mag # unit normal j
      norms[2,ielem]=1.0/mag        # unit normal k

  return nodes, elems, norms



def Wood95Operator(nxin, nyin, dx, dy, cx, cy, lamb, hmax, dim_type, offsetx, offsety, offsetz, noise, halo):
  nx = nxin+(2*halo)
  ny = nyin+(2*halo)
  field = np.zeros((nx,ny))
  dhdx = np.zeros((nx,ny))
  dhdy = np.zeros((nx,ny))
  nodes = np.zeros((3,(nx*ny)))
  elems = np.zeros((4,(nx*ny)), dtype=int)
  norms = np.zeros((3,(nx*ny)))
  in_arr = np.zeros((nx,ny), dtype=int)

  inode = -1
  for i in range(nx):
    x = ((i-halo)*dx + offsetx) - cx
    for j in range(ny):
      inode += 1
      in_arr[i,j]=inode+1 # Fortran indexing
      y = ((j-halo)*dy + offsety) - cy

      if dim_type == 1:
        r = x/lamb
      elif dim_type == 2:
        r = np.sqrt((x/lamb)**2.0 + (y/lamb)**2.0)
      else:
        r = -999.0
      if np.abs(r) <= 0.5:
        field[i,j] = hmax*np.cos(np.pi*r)**2.0
        if noise:
          field[i,j] = field[i,j]+(np.random.rand()*0.1)
      else:
        field[i,j] = 0.0


      # gradients
  for j in range(ny):
    dhdx[:-1,j] = np.diff(field[:,j])/dx
    dhdx[-1,j] = (field[0,j]-field[-1,j])/dx
  for i in range(nx):
    dhdy[i,:-1] = np.diff(field[i,:])/dy
    dhdy[i,-1] = (field[i,0]-field[i,-1])/dy

  inode = -1
  ielem = -1
  for i in range(nx):
    for j in range(ny):
      inode += 1
      ielem += 1
      ip1 = i+1
      jp1 = j+1
      if ip1 > nx-1: ip1 = i
      if jp1 > ny-1: jp1 = j
      nodes[0,inode]=(i-halo)*dx + offsetx
      nodes[1,inode]=(j-halo)*dy + offsety
      nodes[2,inode]=field[i,j] + offsetz
      elems[0,ielem]=in_arr[i,j]     # associated nodes (clockwise)
      elems[1,ielem]=in_arr[ip1,j]
      elems[2,ielem]=in_arr[ip1,jp1]
      elems[3,ielem]=in_arr[i,jp1]
      mag = np.sqrt((dhdx[i,j]*dhdx[i,j]) + (dhdy[i,j]*dhdy[i,j]) + 1.0)
      norms[0,ielem]=-dhdx[i,j]/mag # unit normal i
      norms[1,ielem]=-dhdy[i,j]/mag # unit normal j
      norms[2,ielem]=1.0/mag        # unit normal k

  return nodes, elems, norms


def FlatOperator(hmax, nxin, nyin, dx, dy, offsetx, offsety, halo):
  nx = nxin+(2*halo)
  ny = nyin+(2*halo)
  field = np.zeros((nx,ny))
  nodes = np.zeros((3,(nx*ny)))
  elems = np.zeros((4,(ny*nx)), dtype=int)
  norms = np.zeros((3,(nx*ny)))
  in_arr = np.zeros((nx,ny), dtype=int)

  inode = -1
  for i in range(nx):
    for j in range(ny):
      inode += 1
      in_arr[i,j]=inode+1 # Fortran indexing
      field[i,j] = hmax

  inode = -1
  ielem = -1
  for i in range(nx):
    for j in range(ny):
      inode += 1
      ielem += 1
      ip1 = i+1
      jp1 = j+1
      if ip1 > nx-1: ip1 = i
      if jp1 > ny-1: jp1 = j
      nodes[0,inode]=(i-halo)*dx
      nodes[1,inode]=(j-halo)*dy
      nodes[2,inode]=field[i,j]
      elems[0,ielem]=in_arr[i,j]     # associated nodes (clockwise)
      elems[1,ielem]=in_arr[ip1,j]
      elems[2,ielem]=in_arr[ip1,jp1]
      elems[3,ielem]=in_arr[i,jp1]
      norms[0,ielem]=0.0
      norms[1,ielem]=0.0
      norms[2,ielem]=1.0

  return nodes, elems, norms



def BlocksOperator(cents, hdims, hghts, dx, dy, xsize, ysize, offsetx, offsety, offsetz, halo):

  nx = int(round(xsize/dx))+(2*halo)
  ny = int(round(ysize/dy))+(2*halo)
  dz = min(dx,dy)
  nodes = [[],[],[]]
  elems = [[],[],[],[]]
  norms = [[],[],[]]
  xmax = (nx-1)*dx
  ymax = (ny-1)*dy

  xls = []
  xrs = []
  yls = []
  yrs = []
  ngroup0 = [] # node group0
  ngroups = [] # node groups
  egroup0 = [] # element group0
  egroups = [] # element groups
  nblocks = len(cents)
  inode = -1
  ielem = -1
  nz = []

  for b in range(nblocks):
    ngroups.append([[],[],[],[],[]]) # node index lists for meshgroups (xl, xr, yl, yr, top)
    egroups.append([[],[],[],[],[]]) # elemnt index lists for meshgroups (xl, xr, yl, yr, top)
    xls.append(int(round((cents[b][0] - hdims[b][0]/2.)/dx)))
    xrs.append(int(round((cents[b][0] + hdims[b][0]/2.)/dx)))
    yls.append(int(round((cents[b][1] - hdims[b][1]/2.)/dy)))
    yrs.append(int(round((cents[b][1] + hdims[b][1]/2.)/dy)))
    nz.append(int(hghts[b]/dz)+1)
    if isclose(np.asarray([hghts[b]]), dz*(nz[b]-1)):nz[b] -= 1

  # flat ground meshgroup
  for ii in range(nx):
    for jj in range(ny):
      i=ii-halo
      j=jj-halo
      makenode = True
      makeelem = True
      for b in range(nblocks):
        if (i > xls[b] and i < xrs[b] and j > yls[b] and j < yrs[b]):
          makenode=False
          makeelem=False
      if makenode: inode+=1
      for b in range(nblocks):
          if (i == xls[b] and j >= yls[b] and j < yrs[b]):
            egroups[b][0].append(inode)
            makeelem=False
          if (j == yls[b] and i >= xls[b] and i < xrs[b]):
            egroups[b][2].append(inode)
            makeelem=False
          if (i == xrs[b] and j >= yls[b] and j < yrs[b]):
            egroups[b][1].append(inode)
          if (j == yrs[b] and i >= xls[b] and i < xrs[b]):
            egroups[b][3].append(inode)

          if (i == xls[b] and j >= yls[b] and j <= yrs[b]): ngroups[b][0].append(inode)
          if (i == xrs[b] and j >= yls[b] and j <= yrs[b]): ngroups[b][1].append(inode)
          if (i >= xls[b] and i <= xrs[b] and j == yls[b]): ngroups[b][2].append(inode)
          if (i >= xls[b] and i <= xrs[b] and j == yrs[b]): ngroups[b][3].append(inode)
      if makenode:
        nodes[0].append(i*dx + offsetx)
        nodes[1].append(j*dy + offsety)
        nodes[2].append(0.0)
        ngroup0.append(inode)
        if makeelem:egroup0.append(inode)

  nodemaxx=np.max(nodes[0][:])
  nodeminx=np.min(nodes[0][:])
  nodemaxy=np.max(nodes[1][:])
  nodeminy=np.min(nodes[1][:])
  
  for b in range(nblocks):
    # XL meshgroup
    for k in range(1,nz[b]):
      for j in range(yls[b],yrs[b]+1):
        inode+=1
        nodes[0].append(xls[b]*dx + offsetx)
        nodes[1].append(j*dy + offsety)
        nodes[2].append(k*dz)
        ngroups[b][0].append(inode)
        if j == yls[b]:
          ngroups[b][2].append(inode)
          egroups[b][2].append(inode)
        if j == yrs[b]:
          ngroups[b][3].append(inode)
          egroups[b][3].append(inode)
        if j >= yls[b] and j < yrs[b]:
          egroups[b][0].append(inode)

    for j in range(yls[b],yrs[b]+1):
      inode+=1
      nodes[0].append(xls[b]*dx + offsetx)
      nodes[1].append(j*dy + offsety)
      nodes[2].append(hghts[b])
      ngroups[b][0].append(inode)
      ngroups[b][4].append(inode)
      if j >= yls[b] and j < yrs[b]:
        egroups[b][4].append(inode)
      if j == yls[b]:
        ngroups[b][2].append(inode)
      if j == yrs[b]:
        ngroups[b][3].append(inode)
      
    # XR meshgroup
    for k in range(1,nz[b]):
      for j in range(yls[b],yrs[b]+1):
        inode+=1
        nodes[0].append(xrs[b]*dx + offsetx)
        nodes[1].append(j*dy + offsety)
        nodes[2].append(k*dz)
        ngroups[b][1].append(inode)
        if j == yls[b]:
          ngroups[b][2].append(inode)
        if j == yrs[b]:
          ngroups[b][3].append(inode)
        if j >= yls[b] and j < yrs[b]:
          egroups[b][1].append(inode)

    for j in range(yls[b],yrs[b]+1):
      inode+=1
      nodes[0].append(xrs[b]*dx + offsetx)
      nodes[1].append(j*dy + offsety)
      nodes[2].append(hghts[b])
      ngroups[b][1].append(inode)
      ngroups[b][4].append(inode)
      if j == yls[b]:
        ngroups[b][2].append(inode)
      if j == yrs[b]:
        ngroups[b][3].append(inode)
      

    # YL meshgroup
    for k in range(1,nz[b]):
      for i in range(xls[b]+1,xrs[b]):
        inode+=1
        nodes[0].append(i*dx + offsetx)
        nodes[1].append(yls[b]*dy + offsety)
        nodes[2].append(k*dz)
        ngroups[b][2].append(inode)
        egroups[b][2].append(inode)

    for i in range(xls[b]+1,xrs[b]):
      inode+=1
      nodes[0].append(i*dx + offsetx)
      nodes[1].append(yls[b]*dy + offsety)
      nodes[2].append(hghts[b])
      ngroups[b][2].append(inode)
      ngroups[b][4].append(inode)
      egroups[b][4].append(inode)


    # YR meshgroup
    for k in range(1,nz[b]):
      for i in range(xls[b]+1,xrs[b]):
        inode+=1
        nodes[0].append(i*dx + offsetx)
        nodes[1].append(yrs[b]*dy + offsety)
        nodes[2].append(k*dz)
        ngroups[b][3].append(inode)
        egroups[b][3].append(inode)

    for i in range(xls[b]+1,xrs[b]):
      inode+=1
      nodes[0].append(i*dx + offsetx)
      nodes[1].append(yrs[b]*dy + offsety)
      nodes[2].append(hghts[b])
      ngroups[b][3].append(inode)
      ngroups[b][4].append(inode)


    # TOP meshgroup
    for j in range(yls[b]+1,yrs[b]):
      for i in range(xls[b]+1,xrs[b]):
        inode+=1
        nodes[0].append(i*dx + offsetx)
        nodes[1].append(j*dy + offsety)
        nodes[2].append(hghts[b])
        ngroups[b][4].append(inode)
        egroups[b][4].append(inode)



  # make elements from group lists

  nodes = np.asarray(nodes)
  egroup0 = np.asarray(egroup0)
  egroups = np.asarray(egroups)

# group 0
  ng = np.asarray(ngroup0)
  for igrp in range(len(egroup0)):
    p0 = nodes[0,egroup0[igrp]]
    q0 = nodes[1,egroup0[igrp]]
    p1 = p0 + dx
    q1 = q0 + dy
    if p1>nodemaxx:p1=p0
    if q1>nodemaxy:q1=q0
    idtmp = np.where(isclose(nodes[0,ng],p1))[0]
    idx1 = idtmp[np.where(isclose(nodes[1,ng[idtmp]],q0))][0]
    idx2 = idtmp[np.where(isclose(nodes[1,ng[idtmp]],q1))][0]
    idtmp = np.where(isclose(nodes[0,ng],p0))[0]
    idx3 = idtmp[np.where(isclose(nodes[1,ng[idtmp]],q1))][0]
    elems[0].append(egroup0[igrp]+1)
    elems[1].append(ng[idx1]+1)
    elems[2].append(ng[idx2]+1)
    elems[3].append(ng[idx3]+1)
    norms[0].append(0.0)
    norms[1].append(0.0)
    norms[2].append(1.0)

  for b in range(nblocks):
    # XL
    face=0
    ng = np.asarray(ngroups[b][face])
    for igrp in range(len(egroups[b][face])):
      p0 = nodes[1,egroups[b][face][igrp]]
      q0 = nodes[2,egroups[b][face][igrp]]
      p1 = p0 + dy
      q1 = q0 + dz
      if isclose(np.asarray([p1]),(ysize)):p1=0.0+offsetx
      if q1 > hghts[b]:q1=hghts[b]
      idtmp = np.where(isclose(nodes[1,ng],p1))[0]
      idx1 = idtmp[np.where(isclose(nodes[2,ng[idtmp]],q0))][0]
      idx2 = idtmp[np.where(isclose(nodes[2,ng[idtmp]],q1))][0]
      idtmp = np.where(isclose(nodes[1,ng],p0))[0]
      idx3 = idtmp[np.where(isclose(nodes[2,ng[idtmp]],q1))][0]
      elems[0].append(egroups[b][face][igrp]+1)
      elems[1].append(ng[idx1]+1)
      elems[2].append(ng[idx2]+1)
      elems[3].append(ng[idx3]+1)
      norms[0].append(-1.0)
      norms[1].append(0.0)
      norms[2].append(0.0)
  
    # XR
    face=1
    ng = np.asarray(ngroups[b][face])
    for igrp in range(len(egroups[b][face])):
      p0 = nodes[1,egroups[b][face][igrp]]
      q0 = nodes[2,egroups[b][face][igrp]]
      p1 = p0 + dy
      q1 = q0 + dz
      if isclose(np.asarray([p1]),(ysize)):p1=0.0+offsetx
      if q1 > hghts[b]:q1=hghts[b]
      idtmp = np.where(isclose(nodes[1,ng],p1))[0]
      idx1 = idtmp[np.where(isclose(nodes[2,ng[idtmp]],q0))][0]
      idx2 = idtmp[np.where(isclose(nodes[2,ng[idtmp]],q1))][0]
      idtmp = np.where(isclose(nodes[1,ng],p0))[0]
      idx3 = idtmp[np.where(isclose(nodes[2,ng[idtmp]],q1))][0]
      elems[0].append(egroups[b][face][igrp]+1)
      elems[1].append(ng[idx1]+1)
      elems[2].append(ng[idx2]+1)
      elems[3].append(ng[idx3]+1)
      norms[0].append(1.0)
      norms[1].append(0.0)
      norms[2].append(0.0)
  
    # YL
    face=2
    ng = np.asarray(ngroups[b][face])
    for igrp in range(len(egroups[b][face])):
      p0 = nodes[0,egroups[b][face][igrp]]
      q0 = nodes[2,egroups[b][face][igrp]]
      p1 = p0 + dx
      q1 = q0 + dz
      if isclose(np.asarray([p1]),(xsize)):p1=0.0+offsetx
      if q1 > hghts[b]:q1=hghts[b]
      idtmp = np.where(isclose(nodes[0,ng],p1))[0]
      idx1 = idtmp[np.where(isclose(nodes[2,ng[idtmp]],q0))][0]
      idx2 = idtmp[np.where(isclose(nodes[2,ng[idtmp]],q1))][0]
      idtmp = np.where(isclose(nodes[0,ng],p0))[0]
      idx3 = idtmp[np.where(isclose(nodes[2,ng[idtmp]],q1))][0]
      elems[0].append(egroups[b][face][igrp]+1)
      elems[1].append(ng[idx1]+1)
      elems[2].append(ng[idx2]+1)
      elems[3].append(ng[idx3]+1)
      norms[0].append(0.0)
      norms[1].append(-1.0)
      norms[2].append(0.0)
  
    # YR
    face=3
    ng = np.asarray(ngroups[b][face])
    for igrp in range(len(egroups[b][face])):
      p0 = nodes[0,egroups[b][face][igrp]]
      q0 = nodes[2,egroups[b][face][igrp]]
      p1 = p0 + dx
      q1 = q0 + dz
      if isclose(np.asarray([p1]),(xsize)):p1=0.0+offsetx
      if q1 > hghts[b]:q1=hghts[b]
      idtmp = np.where(isclose(nodes[0,ng],p1))[0]
      idx1 = idtmp[np.where(isclose(nodes[2,ng[idtmp]],q0))][0]
      idx2 = idtmp[np.where(isclose(nodes[2,ng[idtmp]],q1))][0]
      idtmp = np.where(isclose(nodes[0,ng],p0))[0]
      idx3 = idtmp[np.where(isclose(nodes[2,ng[idtmp]],q1))][0]
      elems[0].append(egroups[b][face][igrp]+1)
      elems[1].append(ng[idx1]+1)
      elems[2].append(ng[idx2]+1)
      elems[3].append(ng[idx3]+1)
      norms[0].append(0.0)
      norms[1].append(1.0)
      norms[2].append(0.0)
  
    # TOP
    face=4
    ng = np.asarray(ngroups[b][face])
    for igrp in range(len(egroups[b][face])):
      p0 = nodes[0,egroups[b][face][igrp]]
      q0 = nodes[1,egroups[b][face][igrp]]
      p1 = p0 + dx
      q1 = q0 + dy
      if isclose(np.asarray([p1]),(xsize)):p1=0.0+offsetx
      if isclose(np.asarray([q1]),(ysize)):q1=0.0+offsety
      idtmp = np.where(isclose(nodes[0,ng],p1))[0]
      idx1 = idtmp[np.where(isclose(nodes[1,ng[idtmp]],q0))][0]
      idx2 = idtmp[np.where(isclose(nodes[1,ng[idtmp]],q1))][0]
      idtmp = np.where(isclose(nodes[0,ng],p0))[0]
      idx3 = idtmp[np.where(isclose(nodes[1,ng[idtmp]],q1))][0]
      elems[0].append(egroups[b][face][igrp]+1)
      elems[1].append(ng[idx1]+1)
      elems[2].append(ng[idx2]+1)
      elems[3].append(ng[idx3]+1)
      norms[0].append(0.0)
      norms[1].append(0.0)
      norms[2].append(1.0)
  

  elems = np.asarray(elems)
  norms = np.asarray(norms)
  nodes[2,:]+=offsetz

  return nodes, elems, norms


