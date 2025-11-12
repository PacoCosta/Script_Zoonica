# -*- coding: utf-8 -*-
"""
Created on Thu Nov  6 01:24:39 2025

@author: fcosta
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

############################
# Conservacion de la energía
############################


# Parámetros físicos
g = 9.81     # gravedad (m/s²)
L = 1.0      # longitud (m)
M = 1.0      # masa (kg)
theta0 = 0.3 # ángulo inicial (rad)
omega0 = 0.0 # velocidad angular inicial (rad/s)

# Tiempo
t_max = 10
dt = 0.005
t = np.arange(0, t_max, dt)

# Ecuación diferencial del péndulo
def f1(t, y):
    theta, omega = y
    dtheta = omega
    domega = - (g/L) * np.sin(theta)
    return np.array([dtheta, domega])

# Método de integración RK4
def rk4_step(f1, t, y, dt):
    k1 = f1(t, y)
    k2 = f1(t + dt/2, y + dt*k1/2)
    k3 = f1(t + dt/2, y + dt*k2/2)
    k4 = f1(t + dt, y + dt*k3)
    return y + dt*(k1 + 2*k2 + 2*k3 + k4)/6

# Integración temporal
y = np.zeros((len(t), 2))
y[0] = [theta0, omega0]
for i in range(len(t)-1):
    y[i+1] = rk4_step(f1, t[i], y[i], dt)

theta = y[:,0]
omega = y[:,1]
x = L * np.sin(theta)
y_pos = -L * np.cos(theta)

# Energías
E_k = 0.5 * M * (L * omega)**2
E_p = M * g * L * (1 - np.cos(theta))
E_t = E_k + E_p

# --- ANIMACIÓN ---
fig1, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 8))
plt.subplots_adjust(hspace=0.4)

# Subgráfico 1: Péndulo
ax1.set_xlim(-L*1.2, L*1.2)
ax1.set_ylim(-L*1.2, L*0.2)
ax1.set_aspect('equal')
ax1.set_title("Péndulo simple")
ax1.axis('off')
line, = ax1.plot([], [], 'o-', lw=3, color='royalblue')
trace, = ax1.plot([], [], '-', lw=1, color='gray', alpha=0.5)
trail_x, trail_y = [], []

# Subgráfico 2: Energías
ax2.set_xlim(0, t_max)
ax2.set_ylim(0, 1.2 * max(E_t))
ax2.set_xlabel("Tiempo (s)")
ax2.set_ylabel("Energía (J)")
ax2.set_title("Energías del péndulo simple")
line_Ek, = ax2.plot([], [], color='orange', label='Energía cinética')
line_Ep, = ax2.plot([], [], color='green', label='Energía potencial')
line_Et, = ax2.plot([], [], color='red', label='Energía total', lw=2)
ax2.legend(loc='upper right')

def init1():
    line.set_data([], [])
    trace.set_data([], [])
    line_Ek.set_data([], [])
    line_Ep.set_data([], [])
    line_Et.set_data([], [])
    return line, trace, line_Ek, line_Ep, line_Et

def update(frame):
    # --- PÉNDULO ---
    x0, y0 = 0, 0
    x1, y1 = x[frame], y_pos[frame]
    line.set_data([x0, x1], [y0, y1])
    trail_x.append(x1)
    trail_y.append(y1)
    trace.set_data(trail_x, trail_y)
    # --- ENERGÍA ---
    line_Ek.set_data(t[:frame], E_k[:frame])
    line_Ep.set_data(t[:frame], E_p[:frame])
    line_Et.set_data(t[:frame], E_t[:frame])
    return line, trace, line_Ek, line_Ep, line_Et

ani1 = FuncAnimation(fig1, update, frames=len(t),
                    init_func=init1, interval=dt*1000, blit=True)

plt.show()

#######################################################
#######################################################
#######################################################

# -------------------------
# Dinámica del péndulo
# -------------------------
def f2(t, y):
    theta, omega = y
    dtheta = omega
    domega = -(g/L) * np.sin(theta)
    return np.array([dtheta, domega])

def rk4_step(f2, t, y, dt):
    k1 = f2(t, y)
    k2 = f2(t + dt/2, y + dt*k1/2)
    k3 = f2(t + dt/2, y + dt*k2/2)
    k4 = f2(t + dt, y + dt*k3)
    return y + dt*(k1 + 2*k2 + 2*k3 + k4)/6

# Integración temporal
Y = np.zeros((len(t), 2))
Y[0] = [theta0, omega0]
for i in range(len(t)-1):
    Y[i+1] = rk4_step(f2, t[i], Y[i], dt)

theta = Y[:,0]
omega = Y[:,1]

# -------------------------
# Cálculo de fuerzas
# -------------------------
# Tensión
T = M * (L * omega**2 + g * np.cos(theta))

# Peso (vector)
P_x = np.zeros_like(t)
P_y = -M * g * np.ones_like(t)

# Tensión (vector)
T_x = -T * np.sin(theta)
T_y =  T * np.cos(theta)

# Fuerza neta (T + P)
Fnet_x = T_x + P_x
Fnet_y = T_y + P_y
Fnet_mag = np.sqrt(Fnet_x**2 + Fnet_y**2)

# Posición de la masa
x =  L * np.sin(theta)
y = -L * np.cos(theta)

# -------------------------
# CONFIGURAR FIGURA Y SUBPLOTS
# -------------------------
fig2 = plt.figure(figsize=(7,8))
gs = fig2.add_gridspec(2, 1, height_ratios=[2, 1])
ax11 = fig2.add_subplot(gs[0])  # animación
ax22 = fig2.add_subplot(gs[1])  # gráfico de fuerzas

# --- Subgráfico 1: animación ---
ax11.set_xlim(-1.2*L, 1.2*L)
ax11.set_ylim(-1.2*L, 0.3*L)
ax11.set_aspect('equal')
ax11.axis('off')
ax11.set_title("Péndulo simple — movimiento y fuerzas")

line, = ax11.plot([], [], 'o-', lw=3, color='royalblue')
trace, = ax11.plot([], [], '-', lw=1, color='gray', alpha=0.4)
vec_T, = ax11.plot([], [], color='purple', lw=2, label='Tensión T')
vec_P, = ax11.plot([], [], color='red', lw=2, label='Peso mg')
vec_F, = ax11.plot([], [], color='green', lw=2, label='Fuerza neta')
ax11.legend(loc='upper right')

trail_x, trail_y = [], []

# --- Subgráfico 2: gráfico de fuerzas ---
ax22.set_xlim(0, t_max)
ax22.set_ylim(0, 1.2 * max(T))
ax22.set_xlabel("Tiempo (s)")
ax22.set_ylabel("Fuerza (N)")
ax22.set_title("Magnitud de las fuerzas vs tiempo")
line_T, = ax22.plot([], [], color='purple', lw=2, label='Tensión')
line_P, = ax22.plot([], [], color='red', lw=2, label='Peso')
line_F, = ax22.plot([], [], color='green', lw=2, label='F_{net}')
# ax22.legend(loc='upper right')

# -------------------------
# FUNCIONES DE ANIMACIÓN
# -------------------------
def init2():
    for l in [line, trace, vec_T, vec_P, vec_F, line_T, line_P, line_F]:
        l.set_data([], [])
    return line, trace, vec_T, vec_P, vec_F, line_T, line_P, line_F

def update2(frame):
    # Posición del péndulo
    x0, y0 = 0, 0
    x1, y1 = x[frame], y[frame]
    line.set_data([x0, x1], [y0, y1])

    # Trayectoria
    trail_x.append(x1)
    trail_y.append(y1)
    trace.set_data(trail_x[-200:], trail_y[-200:])  # dejar solo los últimos puntos

    # Escala de vectores
    scale = 0.25
    vec_T.set_data([x1, x1 + T_x[frame]*scale/T.max()],
                   [y1, y1 + T_y[frame]*scale/T.max()])
    vec_P.set_data([x1, x1 + P_x[frame]*scale/g],
                   [y1, y1 + P_y[frame]*scale/g])
    vec_F.set_data([x1, x1 + Fnet_x[frame]*scale/T.max()],
                   [y1, y1 + Fnet_y[frame]*scale/T.max()])

    # Actualizar gráfico inferior
    line_T.set_data(t[:frame], T[:frame])
    line_P.set_data(t[:frame], np.ones(frame)*M*g)
    line_F.set_data(t[:frame], Fnet_mag[:frame])

    return line, trace, vec_T, vec_P, vec_F, line_T, line_P, line_F

ani2 = FuncAnimation(fig2, update2, frames=len(t),
                    init_func=init2, interval=dt*1000, blit=True)

plt.tight_layout()
plt.show()
