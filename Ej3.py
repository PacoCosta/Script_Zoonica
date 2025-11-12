# -*- coding: utf-8 -*-
"""
Created on Wed Nov 12 10:15:13 2025

@author: fcosta
"""

import numpy as np
import matplotlib.pyplot as plt

# ============================
# Input parameters
# ============================
theta = np.radians(np.linspace(0,360,361))
alpha = np.radians(np.linspace(0, 360, 361))   


# Grid
alpha_mesh, theta_mesh = np.meshgrid(alpha, theta)

# ============================
# Implementation Eq.(8)
# ============================
Te = (np.sin(alpha_mesh) * np.cos(theta_mesh) - np.cos(theta_mesh) * np.cos(alpha_mesh)**2)
# Te = np.cos(theta_mesh) * (np.sin(alpha_mesh) - np.cos(alpha_mesh)**2)

# ============================
# Figure 5
# ============================
fig1 = plt.figure(figsize=(8,6))
ax = fig1.add_subplot(121, projection='3d')
surf = ax.plot_surface(np.degrees(alpha_mesh), np.degrees(theta_mesh), Te, cmap='plasma')
ax.set_xlabel(r'$\alpha$ [°]',fontsize=20)
ax.set_ylabel(r'$\theta$ [°]',fontsize=20)
ax.set_zlabel('T [-]',fontsize=20)
ax.set_title('Dot Product',fontsize=20)

ax = fig1.add_subplot(122)
im = ax.contourf(Te, extent=[ 0, 359, 0, 359], origin='lower', cmap='plasma')
ax.set_xlabel(r'$\alpha$ [°]')
ax.set_ylabel(r'$\theta$ [°]')
ax.set_title(r'M/m( $\alpha$,$\theta$)')
fig1.colorbar(im, ax=ax, orientation='vertical', shrink=0.8)

plt.show()


# ============================
# Figure 6
# ============================
fig2, ax2 = plt.subplots(1, 2, figsize=(10, 8))


# (c) 
theta_idx = np.argmin((180))
ax2[0].plot(np.degrees(alpha), Te[179])
ax2[0].set_title(r'Corte en $\theta$ = 180°',fontsize=20)
ax2[0].set_xlabel(r'$\alpha$ [°]',fontsize=20)
ax2[0].set_ylabel('T [-]',fontsize=20)
ax2[0].grid(True)

# (d) 
alpha_idx = np.argmin(36.37)
ax2[1].plot(np.degrees(theta), Te[:, alpha_idx])
ax2[1].set_title(r'Corte en $\alpha$ = 36.37°',fontsize=20)
ax2[1].set_xlabel(r'$\theta$ [°]',fontsize=20)
ax2[1].set_ylabel('T [-]',fontsize=20)
ax2[1].grid(True)

plt.tight_layout()
plt.show()
