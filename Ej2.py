# -*- coding: utf-8 -*-
"""
Created on Wed Nov 12 01:02:48 2025

@author: fcosta
"""
import numpy as np
import matplotlib.pyplot as plt

# ============================
# Input parameters
# ============================
g     = 9.81                # gravedad [m/s²]
m     = 1.0                 # masa del péndulo [kg]
M     = 1.0
r1    = 1.0
r2    = 1.0
alpha = np.radians(np.linspace(0,90,100))


# ============================
# Figure 3
# ============================
Fd       = m*g*(2*np.sin(alpha)+1)*np.cos(alpha)
dfdalpha = m*g*((2*np.cos(alpha)**2)-(2*np.sin(alpha)**2)-np.sin(alpha))

fig1,ax1 = plt.subplots()
ax1.plot(np.degrees(alpha),Fd)
ax1.plot(np.degrees(alpha),dfdalpha)
ax1.set_ylabel(r'$F_d$',fontsize=20)
ax1.set_xlabel(r'$\alpha$ [°]',fontsize=20)
ax1.grid(axis='both')
ax1.tick_params(axis='both', labelsize=20)

# ============================
# Figure 4
# ============================
mu = 0.6
Fd2       = m*g*(2*np.sin(alpha)+1)*np.cos(alpha)-mu*g*(m+M)
dfdalpha2 = m*g*((2*np.cos(alpha)**2)-(2*np.sin(alpha)**2)-np.sin(alpha))

fig2,ax2 = plt.subplots()
ax2.plot(np.degrees(alpha),Fd2)
ax2.plot(np.degrees(alpha),dfdalpha2)
ax2.set_ylabel(r'$F_d$',fontsize=20)
ax2.set_xlabel(r'$\alpha$ [°]',fontsize=20)
ax2.grid(axis='both')
ax2.tick_params(axis='both', labelsize=20)

ax2.text(75,15, '$\mu$ = {}'.format(mu), fontsize = 16, bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.4) ) 

 


# ============================
# Figure 6
# ============================

# 1. Inputs

theta = np.radians(np.arange(0, 400, 1))
alpha = np.radians(np.arange(0, 400, 1))

#================================
# 2. M/m matrix simplifying mats as follows:

# e^(j theta) + e^(-j theta)       =  2 cos(theta)
# (e^(-j alpha) - e^(j alpha))^2   = (-2 j sin^2(theta)) = - 4 sin^2(alpha)

# Therefore:
# Mm = 1 + 4 cos(theta) sin^2(alpha)
#================================

Mm = np.zeros((len(alpha), len(theta)))

for i,a in enumerate(alpha):
    for j,t  in enumerate(theta):
        Mm[i,j] = 1 + 4 * np.cos(t) * np.sin(a)**2

#================================
# 3. Figure 6 (1)
#================================
A,T  = np.meshgrid(np.degrees(theta), np.degrees(alpha))
fig  = plt.figure(figsize=(10, 6))
ax   = fig.add_subplot(121, projection='3d')
surf = ax.plot_surface(T, A, Mm, cmap='plasma')

ax.set_xlabel(r'$\alpha$ [°]',fontsize=20)
ax.set_ylabel(r'$\theta$ [°]',fontsize=20)
ax.set_zlabel('M/m',fontsize=20)
ax.set_title(r' M/m($\alpha$, $\theta$)',fontsize=20)

#================================
# Figure 6 (2)
#================================
ax2 = fig.add_subplot(122)
im  = ax2.imshow(Mm, extent=[0, 400, 0, 400],
                origin='lower', cmap='plasma', aspect='auto')
ax2.set_xlabel(r'$\theta$ [°]',fontsize=20)
ax2.set_ylabel(r'$\alpha$ [°]',fontsize=20)
# ax2.set_title('Mapa de colores M/m(θ, α)')
fig.colorbar(im, ax=ax2, orientation='vertical', shrink=0.8)

plt.tight_layout()
plt.show()

