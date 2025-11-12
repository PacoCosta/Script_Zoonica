# -*- coding: utf-8 -*-
"""
Created on Tue Nov 11 09:25:29 2025

@author: fcosta
"""
import numpy as np
import matplotlib.pyplot as plt


# ============================
# Input parameters
# ============================
g     = 9.81                
m     = 1.0                 # masa del péndulo [kg]
M     = 1.0                 # masa colgante [kg]
r1    = 1.0                 # distancia fulcro -> masa m [m]
r2    = 1.0                 # distancia fulcro -> masa M [m]
l     = 1                   # longitud del péndulo [m]
alpha = 1000/(3600/(2*np.pi))   # angular acceleration 
 

# ============================
# Figure 3
# ============================

theta = np.radians(np.linspace(0,180,1000))
w     = np.sqrt((2*g/l)*(1-np.cos(theta)))

fig1,ax0 = plt.subplots()
ax0.plot(np.degrees(theta),w)
ax0.set_ylabel(r'$\omega$ [rad/s]',fontsize=20)
ax0.set_xlabel(r'$\theta$ [°]',fontsize=20)
ax0.grid(axis='both')
ax0.tick_params(axis='both', labelsize=20)


# ============================
#Figure 6
# ============================

wf = (np.linspace(0,5000,5000))
w0 = 0
dE = (wf**2-w0**2)*m*l**2/2

fig2,ax2 = plt.subplots()
ax2.plot(wf,dE)
ax2.set_ylabel(r'$\Delta$E [J]',fontsize=20)
ax2.set_xlabel(r'$\omega$ [rad/s]',fontsize=20)
ax2.grid(axis='both')
ax2.tick_params(axis='both', labelsize=20)



# ============================
# Figure 8
# ============================
theta0=np.radians([30,90,180])
theta=np.radians(np.linspace(-180,180,1000))

wf=(np.sqrt((2*g/l)*(1-np.cos(theta0))))

fig3,ax3 = plt.subplots()
Mm=[]
for i in (range(len(wf))):
    Mm=((r1/r2)*((l*(wf[i]**2)*np.cos(theta)/g)+1))
    ax3.plot(np.degrees(theta),Mm,label=r'$\theta_0$ ={}°'.format(np.round(np.degrees(theta0[i]),2)))

ax3.set_ylabel('M/m',fontsize=20)
ax3.set_xlabel(r'$\theta$ [°]',fontsize=20)
ax3.grid(axis='both')
ax3.tick_params(axis='both', labelsize=20)
plt.legend()
plt.show


# ============================
# Figure 9
# ============================

theta=np.radians(np.linspace(-180,180,1000))
theta0 = np.radians(30)
# w0=(1/(2*np.pi))*np.sqrt(g/l)

w0=np.sqrt((2*g/l)*(1-np.cos(theta0)))

wf=([2*np.pi*100/60,2*np.pi*150/60,2*np.pi*200/60])

fig3,ax3 = plt.subplots()
Mm=[]
for i in (range(len(wf))):
    Mm=((r1/r2)*((l*(wf[i]**2)*np.cos(theta)/g)+1))
    ax3.plot(np.degrees(theta),Mm,label=r'$\theta_0$ ={} [rpm]'.format(np.round((wf[i]*60/((2*np.pi))),2)))

ax3.set_ylabel('M/m',fontsize=20)
ax3.set_xlabel(r'$\theta$ [°]',fontsize=20)
ax3.grid(axis='both')
ax3.tick_params(axis='both', labelsize=20)
plt.legend()
plt.show

# ============================
# Figure 10
# ============================

theta=np.radians(np.linspace(-180,180,1000))
w0=(1/(2*np.pi))*np.sqrt(g/l)
wf=([2*np.pi*1000/60,2*np.pi*2000/60])

fig4,ax4 = plt.subplots()
Mm=[]
for i in (range(len(wf))):
    Mm=((r1/r2)*((l*(wf[i]**2)*np.cos(theta)/g)+1))
    ax4.plot(np.degrees(theta),Mm,label=r'$\theta_0$ ={} [rpm]'.format(np.round((wf[i]*60/((2*np.pi))),2)))

ax4.set_ylabel('M/m',fontsize=20)
ax4.set_xlabel(r'$\theta$ [°]',fontsize=20)
ax4.grid(axis='both')
ax4.tick_params(axis='both', labelsize=20)
plt.legend()
plt.show


# ============================
# Figure 11
# ============================

theta = np.radians(np.linspace(0,360,1000))
# w0    = (1/(2*np.pi))*np.sqrt(g/l)
w0    = 60/2*np.pi # 2*np.pi*1000/60 #



Mm = (r1/r2) * (l*(w0**2)*np.cos(theta) + 2*l*theta*alpha*np.cos(theta) + g)
    

fig5,ax5 = plt.subplots()
ax5.plot(np.degrees(theta),Mm)
ax5.set_ylabel('M/m',fontsize=20)
ax5.set_xlabel(r'$\theta$ [°]',fontsize=20)
ax5.grid(axis='both')
ax5.tick_params(axis='both', labelsize=20)
# plt.legend()
plt.show


# ============================
#Figure 12
# ============================
alpha1=[100/(3600/(2*np.pi)),1000/(3600/(2*np.pi)),2000/(3600/(2*np.pi))]
theta0 = np.radians(-180)

theta = np.radians(np.linspace(0,750,1000))
wf    = (2*g/l)*(1-np.cos(theta0))


fig5,ax5 = plt.subplots()
for i in range(len(alpha1)):
    w0    =  (-2*theta*alpha1[i] + wf)

    Fn = -((r1/r2) * (l* w0 + 2*l*theta*alpha1[0] - 4*g) * np.cos(theta) * m)
    Ft = g*np.sin(theta) + alpha1[i]*l * m



    ax5.plot(np.degrees(theta),Fn,label=r'$F_n$ - {}rpm'.format(alpha1[i]*(3600/(2*np.pi))))
    ax5.plot(np.degrees(theta),Ft,label=r'$F_t$  - {}rpm'.format(alpha1[i]*(3600/(2*np.pi))))

ax5.set_ylabel('F [N]',fontsize=20)
ax5.set_xlabel(r'$\theta$ [°]',fontsize=20)
ax5.grid(axis='both')
ax5.tick_params(axis='both', labelsize=20)
plt.legend()
plt.show



# ============================
# Figure 14
# ============================

theta0 = np.radians(-180)

theta = np.radians(np.linspace(0,750,1000))
wf    = (2*g/l)*(1-np.cos(theta0))

M=((r1/r2)*((l*(wf**2)*np.cos(theta)/g)+1))
fig5,ax5 = plt.subplots()

w0    =  (-2*theta*alpha + wf)

Fn1 = -((r1/r2) * ((l*(w0) - 4*g) * np.cos(theta) + alpha) * m)
Fn2 = -((r1/r2) * ((l*(w0) - 4*g) * np.cos(theta) ) * M)

# Ft = (g*np.sin(theta) + alpha)*2*m



ax5.plot(np.degrees(theta),Fn1,label=r'$F_1$ - {}rpm'.format(alpha*(3600/(2*np.pi))))
# ax5.plot(np.degrees(theta),Fn2,label=r'$F_2$  - {}rpm'.format(alpha*(3600/(2*np.pi))))

ax5.set_ylabel('F [N]',fontsize=20)
ax5.set_xlabel(r'$\theta$ [°]',fontsize=20)
ax5.grid(axis='both')
ax5.tick_params(axis='both', labelsize=20)
plt.legend()
plt.show



