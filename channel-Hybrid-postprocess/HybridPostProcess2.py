#!/usr/bin/env python
# This file reads the binary restart file
import os
import sys
import numpy as np
#from matplotlib import pyplot
import matplotlib.pyplot as plt
import pdb
import os.path
import  scipy.io
import math
from scipy.fftpack import fft,ifft

Version="0.1"
gama=1.4
R=0.31745

input_binary_folder='.'
#filename = 'casename_n0025743_t9.000e+01.res'
#filename='CG_Re1000_64_n002858_t1.000e+01.res'
filename='casename_n0025743_t9.000e+01.res'

# Header with integer type
data_type_string = "(12)i4, "
file_header_int = np.fromfile(input_binary_folder + '/'+ filename, dtype = np.dtype(data_type_string))[0]

precision_bytes, version, time_step, Ni, Nj, Nk, nvars = file_header_int[:7]
Ntot = Ni * Nj * Nk
print file_header_int[:7]
#precision_bytes, version, time_step, Ni,    Nj, Nk, nvars
#[    4           50001      25743   512   240   256     5]

# Header with double type
data_type_string += "(1)f8, ("+repr(Ni)+")f8, ("+repr(Nj)+")f8, ("+repr(Nk)+")f8, "
file_header_doub = np.fromfile(input_binary_folder + '/'+ filename, dtype = np.dtype(data_type_string))[0]
time = file_header_doub[1]
x = file_header_doub[2]
y = file_header_doub[3]
z = file_header_doub[4]
#print 'time=',time,'x=',x,'y=',y,'z=',z
#x=0, 18.80320312
#y=-0.99893432,0.99893432
#z=0,6.25546875
print "x-begin=", x[0],"x-end",x[Ni-1]
print "y-begin=", y[0],"y-end",y[Nj-1]
print "z-begin=", z[0],"z-end",z[Nk-1]

# Properties
data_type_string += "("+repr(Ntot * nvars)+")f"+repr(precision_bytes)+", "

print "Reading files..."
prop = np.fromfile(input_binary_folder + '/'+ filename, dtype = np.dtype(data_type_string))[0][5]

rho = np.reshape(prop[:Ntot], [Ni, Nj, Nk])
rho_u   = np.reshape(prop[Ntot:Ntot*2], [Ni, Nj, Nk])
rho_v   = np.reshape(prop[Ntot*2:Ntot*3], [Ni, Nj, Nk])
rho_w   = np.reshape(prop[Ntot*3:Ntot*4], [Ni, Nj, Nk])
rho_E   = np.reshape(prop[Ntot*4:Ntot*5], [Ni, Nj, Nk])
print "rho-begin=",rho[0,0,0],"rho-end=",rho[Ni-1,Nj-1,Nk-1]
print "rho-u-begin=",rho_u[0,0,0],"rho-u-end=",rho_u[Ni-1,Nj-1,Nk-1]
print "rho-v-begin=",rho_v[0,0,0],"rho-v-end=",rho_v[Ni-1,Nj-1,Nk-1]
print "rho-w-begin=",rho_w[0,0,0],"rho-w-end=",rho_w[Ni-1,Nj-1,Nk-1]
print "rho-e-begin=",rho_E[0,0,0],"rho-e-end=",rho_E[Ni-1,Nj-1,Nk-1]
print "************************"
#print "test",rho[Ni,Nj,Nk]
print "************************"
#output the file in matlab compatible format
#scipy.io.savemat('./casename_n0025743_t9.000e+01.mat', mdict={'rho':rho, "rhoU":rho_u, "rhoV":rho_v,"rhoW":rho_w,"rhoE":rho_E})

print "finish reading files .."
u=rho_u/rho
v=rho_v/rho
w=rho_w/rho
#problem 1
u_aver=np.average(np.average(u,axis=0),axis=1)
v_aver=np.average(np.average(v,axis=0),axis=1)
w_aver=np.average(np.average(w,axis=0),axis=1)
fig, ax = plt.subplots()
ax.plot(y, u_aver, 'r', label='U mean')
ax.plot(y, v_aver, 'g:', label='V mean')
ax.plot(y, w_aver, 'k',  label='W maen')
legend = ax.legend(loc='center', shadow=True, fontsize='x-large')

plt.show()
#problem 2  Reynolds stress uu,vv,ww,uv,uw,vw
uu=np.zeros(Nj)
vv=np.zeros(Nj)
ww=np.zeros(Nj)
uv=np.zeros(Nj)
uw=np.zeros(Nj)
vw=np.zeros(Nj)
for j in range(0,Nj):
    uu[j]=np.average(pow((u[:,j,:]-u_aver[j]),2))
    vv[j]=np.average(pow((v[:,j,:]-v_aver[j]),2))
    ww[j]=np.average(pow((w[:,j,:]-w_aver[j]),2))
    uv[j]=np.average((u[:,j,:]-u_aver[j])*(v[:,j,:]-v_aver[j]))
    uw[j]=np.average((u[:,j,:]-u_aver[j])*(w[:,j,:]-w_aver[j]))
    vw[j]=np.average((v[:,j,:]-v_aver[j])*(w[:,j,:]-w_aver[j]))


fig, ax = plt.subplots()
ax.plot(y, uu, 'r', label='uu')
ax.plot(y, vv, 'g:', label='vv')
ax.plot(y, ww, 'k',  label='ww')
ax.plot(y, uv, 'ro', label='uv')
ax.plot(y, uw, 'go:', label='uw')
ax.plot(y, vw, 'ko',  label='vw')
legend = ax.legend(loc='upper center', shadow=True, fontsize='x-large')

plt.show()        
#problem 3
c_2=(gama*(gama-1)*( rho_E/rho-0.5*( pow(u,2)+pow(v,2)+pow(w,2) ) ))
c=pow(c_2,0.5)

c_aver=np.average(np.average(c,axis=0),axis=1)
Ma_aver=u_aver/c_aver
fig, ax = plt.subplots()
ax.plot(y, Ma_aver, 'r', label='Ma')
legend = ax.legend(loc='lower center', shadow=True, fontsize='x-large')
plt.show() 
#problem 4
e=rho_E/rho-0.5*( pow(u,2)+pow(v,2)+pow(w,2) )
p=(gama-1)*rho*e
T=pow(p/rho/R, 0.5)
T_aver=np.average(np.average(T,axis=0),axis=1)

fig, ax = plt.subplots()
ax.plot(y, T_aver, 'r', label='T')
legend = ax.legend(loc='lower center', shadow=True, fontsize='x-large')
plt.show() 



#problem 5 Compute the average temperature gradient at the bottom and top wall.
#************************
grad_T_u = np.zeros([Ni,Nk]) 
grad_T_d = np.zeros([Ni,Nk]) 
grad_T_d[:,:]=(T[:,1,:]-T[:,0,:])/(y[1]-y[0])
grad_T_u[:,:]=(T[:,Nj-1,:]-T[:,Nj-2,:])/(y[Nj-1]-y[Nj-2])
grad_T_aver_d=np.sum(grad_T_d)/(Ni*Nk)
grad_T_aver_u=np.sum(grad_T_u)/(Ni*Nk)
print "T gradient at down wall=", grad_T_aver_d
print "T gradient at up wall =",  grad_T_aver_u

#6. Compute the average centerline and bulk Mach number.
print "mean centerline Ma number=", Ma_aver[Nj/2]

u_m=np.trapz(u_aver,y)/2.0
a_w=pow(gama*R*0.5*(T_aver[0]+T_aver[Nj-1]), 0.5)
print "u mean=",u_m,"a_w=",a_w
print "Ma_b=",u_m/a_w

#problem 7 
#Compute the planar average dilatation and show, at each y location, the minimum and maximum dilatation.
dudx=np.zeros([Ni,Nj,Nk])
dvdy=np.zeros([Ni,Nj,Nk])
dwdz=np.zeros([Ni,Nj,Nk])

dx=x[1]-x[0]
dz=z[1]-z[0]
dudx[1:-1,:,:]=(u[2:,:,:]-u[:-2,:,:])*0.5/dx
dwdz[:,:,1:-1]=(w[:,:,2:]-w[:,:,:-2])*0.5/dz
for j in range(Nj-1):
    dy = y[j+1]-y[j]
    dvdy[:,j,:]=(v[:,j+1,:]-v[:,j,:])/dy
dvdy[:,Nj-1,:]=dvdy[:,Nj-2]

dudx[0,:,:]=dudx[1,:,:]
dudx[Ni-1,:,:]=dudx[Ni-2,:,:]
dvdy[:,0,:]=dvdy[:,1,:]
dvdy[:,Nj-1,:]=dvdy[:,Nj-2,:]
dwdz[:,:,0]=dwdz[:,:,1]
dwdz[:,:,Nk-1]=dwdz[:,:,Nk-2]

dilatation=dudx+dvdy+dwdz
dilatation_aver=np.zeros(Nj)
for j in range(0,Nj):
    dilatation_aver[j]=np.sum(dilatation[:,j,:])
dilatation_aver=dilatation_aver/float(Ni*Nk)

fig, ax = plt.subplots()
ax.plot(y, dilatation_aver, 'r', label='dilatation_aver')
legend = ax.legend(loc='lower center', shadow=True, fontsize='x-large')
plt.show() 

temp1=dilatation.max(axis=0)
temp2=temp1.max(axis=1)
fig, ax = plt.subplots()
ax.plot(y, temp2, 'r', label='dilatation_max')
legend = ax.legend(loc='upper center', shadow=True, fontsize='x-large')
plt.show() 

temp1=dilatation.min(axis=0)
temp2=temp1.min(axis=1)
fig, ax = plt.subplots()
ax.plot(y, temp2, 'r', label='dilatation_min')
legend = ax.legend(loc='lower center', shadow=True, fontsize='x-large')
plt.show() 


#problem 8
#8. Plot the spectrum of the turbulent kinetic energy and dilatation 
#at given y location.
#Ekx_kinetic energy
temp=pow(u[:,Nj/2,Nk/2],2)+pow(v[:,Nj/2,Nk/2],2)+pow(w[:,Nj/2,Nk/2],2)
yy=fft(temp)                     #fft
yf1=abs(yy)/Ni           #uniform
yf2 = yf1[range(int(Ni/2))]  #half the domain due to symmetry

xf  = np.arange(Ni)        # wave number
xf2 = xf[range(int(Ni/2))]/3  #half of the wave number

#plt.subplot(231)
plt.plot(xf2,yf2,'b')
plt.title('Streamwise spectrum of the turbulent kinetic energy)',fontsize=10,color='#F08080')
plt.yscale('log')
plt.xscale('log')
legend = ax.legend(loc='lower left', shadow=True, fontsize='x-large')
plt.show()


#Ekx_u

temp=u[:,Nj/2,Nk/2]
yy=fft(temp)                     #fft
yf1=abs(yy)/Ni           #uniform
yf2 = yf1[range(int(Ni/2))]  #half the domain due to symmetry

xf  = np.arange(Ni)        # wave number
xf2 = xf[range(int(Ni/2))]/3  #half of the wave number

#plt.subplot(231)
plt.plot(xf2,yf2,'b')
plt.title('Ekx_uu)',fontsize=10,color='#F08080')
plt.yscale('log')
plt.xscale('log')
legend = ax.legend(loc='lower left', shadow=True, fontsize='x-large')
plt.show()
#Ekx_v

temp=v[:,Nj/2,Nk/2]
yy=fft(temp)                     #fft
yf1=abs(yy)/Ni           #uniform
yf2 = yf1[range(int(Ni/2))]  #half the domain due to symmetry

xf  = np.arange(Ni)        # wave number
xf2 = xf[range(int(Ni/2))]/3  #half of the wave number

#plt.subplot(231)
plt.plot(xf2,yf2,'b')
plt.title('Ekx_vv)',fontsize=10,color='#F08080')
plt.yscale('log')
plt.xscale('log')
legend = ax.legend(loc='lower left', shadow=True, fontsize='x-large')
plt.show()
#Ekx_w

temp=w[:,Nj/2,Nk/2]
yy=fft(temp)                     #fft
yf1=abs(yy)/Ni           #uniform
yf2 = yf1[range(int(Ni/2))]  #half the domain due to symmetry

xf  = np.arange(Ni)        # wave number
xf2 = xf[range(int(Ni/2))]/3  #half of the wave number

#plt.subplot(231)
plt.plot(xf2,yf2,'b')
plt.title('Ekx_ww)',fontsize=10,color='#F08080')
plt.yscale('log')
plt.xscale('log')
legend = ax.legend(loc='lower left', shadow=True, fontsize='x-large')
plt.show()
#Ekx_dilatation

temp=dilatation[:,Nj/2,Nk/2]
yy=fft(temp)                     #fft
yf1=abs(yy)/Ni           #uniform
yf2 = yf1[range(int(Ni/2))]  #half the domain due to symmetry

xf  = np.arange(Ni)        # wave number
xf2 = xf[range(int(Ni/2))]/3  #half of the wave number

#plt.subplot(231)
#fig=plt.figure(figsize=(3,3))
plt.plot(xf2,yf2,'b')
plt.title('Streamwise spectrum of dilatation)',fontsize=10,color='#F08080')
plt.yscale('log')
plt.xscale('log')
legend = ax.legend(loc='lower left', shadow=True, fontsize='x-large')
plt.show()
#fig.savefig('1.png',dpi=300)

#problem 9 compute the grid spacing to Kolmogorov length scale
#center line Ni/2,Nj/2,Nk/2
nu=0.00005
dy=y[Nj/2]-y[Nj/2-1]
dudx_11=(u[Ni/2+1,Nj/2,Nk/2]-u[Ni/2-1,Nj/2,Nk/2])/2/dx
dudx_12=(u[Ni/2,Nj/2+1,Nk/2]-u[Ni/2,Nj/2-1,Nk/2])/2/dy
dudx_13=(u[Ni/2,Nj/2,Nk/2+1]-u[Ni/2,Nj/2,Nk/2-1])/2/dz

dudx_21=(v[Ni/2+1,Nj/2,Nk/2]-v[Ni/2-1,Nj/2,Nk/2])/2/dx
dudx_22=(v[Ni/2,Nj/2+1,Nk/2]-v[Ni/2,Nj/2-1,Nk/2])/2/dy
dudx_23=(v[Ni/2,Nj/2,Nk/2+1]-v[Ni/2,Nj/2,Nk/2-1])/2/dz

dudx_31=(w[Ni/2+1,Nj/2,Nk/2]-w[Ni/2-1,Nj/2,Nk/2])/2/dx
dudx_32=(w[Ni/2,Nj/2+1,Nk/2]-w[Ni/2,Nj/2-1,Nk/2])/2/dy
dudx_33=(w[Ni/2,Nj/2,Nk/2+1]-w[Ni/2,Nj/2,Nk/2-1])/2/dz
dissipation=pow(dudx_11,2)+pow(dudx_12,2)+pow(dudx_13,2) +\
            pow(dudx_21,2)+pow(dudx_22,2)+pow(dudx_23,2) +\
            pow(dudx_31,2)+pow(dudx_32,2)+pow(dudx_33,2)
kolmogorov_scale=pow(pow(nu,2)/dissipation, 0.25)
resolution_x=dx/kolmogorov_scale
resolution_y=dy/kolmogorov_scale
resolution_z=dz/kolmogorov_scale
print "resolution in x,y,z"
print resolution_x,resolution_y,resolution_z

















