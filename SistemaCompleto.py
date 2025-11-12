import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from math import sin, cos

# ============================
# PARÁMETROS DEL SISTEMA
# ============================
g = 9.81       # gravedad [m/s²]
m = 1.0        # masa del péndulo [kg]
M = 1.0        # masa colgante [kg]
r1 = 1.0       # distancia fulcro -> péndulo [m]
r2 = 1.0       # distancia fulcro -> masa M [m]
l = 0.5        # longitud del péndulo [m]

# ============================
# CONDICIONES INICIALES
# ============================
phi0 = 0.0
phi_dot0 = 0.0
theta0 = 0.1      # ángulo inicial pequeño [rad]
theta_dot0 = 0.0

# ============================
# TIEMPO DE SIMULACIÓN
# ============================
t_max = 20
dt = 0.005
t = np.arange(0, t_max, dt)
n = len(t)

# ============================
# ECUACIONES DINÁMICAS NO LINEALES
# ============================
def accelerations(phi, phidot, theta, thetadot):
    A11 = m * l**2
    A12 = m * r1 * l * np.sin(phi - theta)
    A21 = A12
    A22 = m * r1**2 + M * r2**2

    rhs1 = - m * r1 * l * np.cos(phi - theta) * phidot**2 - m * g * l * np.sin(theta)
    rhs2 = m * r1 * l * np.cos(phi - theta) * thetadot**2 - (M * r2 - m * r1) * g * np.cos(phi)

    A = np.array([[A11, A12],[A21, A22]])
    rhs = np.array([rhs1, rhs2])
    sol = np.linalg.solve(A, rhs)
    ddtheta, ddphi = sol[0], sol[1]
    return ddphi, ddtheta

def derivs(y):
    phi, phidot, theta, thetadot = y
    ddphi, ddtheta = accelerations(phi, phidot, theta, thetadot)
    return np.array([phidot, ddphi, thetadot, ddtheta])

# ============================
# INTEGRACIÓN (RK4)
# ============================
state = np.zeros((n, 4))
state[0] = np.array([phi0, phi_dot0, theta0, theta_dot0])

for i in range(n-1):
    y = state[i]
    k1 = derivs(y)
    k2 = derivs(y + 0.5*dt*k1)
    k3 = derivs(y + 0.5*dt*k2)
    k4 = derivs(y + dt*k3)
    state[i+1] = y + (dt/6.0)*(k1 + 2*k2 + 2*k3 + k4)

phi, phi_dot, theta, theta_dot = state[:,0], state[:,1], state[:,2], state[:,3]

# ============================
# POSICIONES Y VELOCIDADES
# ============================
xP = -r1 * np.cos(phi)
yP = -r1 * np.sin(phi)
xB = xP + l * np.sin(theta)
yB = yP - l * np.cos(theta)
xA =  r2 * np.cos(phi)
yA =  r2 * np.sin(phi)
xM = xA
yM = yA - 0.5

# Velocidades
xP_dot = r1 * np.sin(phi) * phi_dot
yP_dot = -r1 * np.cos(phi) * phi_dot
xB_dot = xP_dot + l * np.cos(theta) * theta_dot
yB_dot = yP_dot + l * np.sin(theta) * theta_dot
xA_dot = -r2 * np.sin(phi) * phi_dot
yA_dot =  r2 * np.cos(phi) * phi_dot
xM_dot = xA_dot
yM_dot = yA_dot

# ============================
# ENERGÍAS
# ============================
E_kin = 0.5*m*(xB_dot**2 + yB_dot**2) + 0.5*M*(xM_dot**2 + yM_dot**2)
E_pot = m*g*yB + M*g*yM
E_total = E_kin + E_pot

# ============================
# TENSIONES Y TORQUES
# ============================
T_bob = np.zeros(n)
T_M = np.zeros(n)
tau_bob = np.zeros(n)
tau_M = np.zeros(n)

for i in range(n):
    phi_i, phidot_i, theta_i, thetadot_i = phi[i], phi_dot[i], theta[i], theta_dot[i]
    ddphi, ddtheta = accelerations(phi_i, phidot_i, theta_i, thetadot_i)
    
    # Aceleraciones
    xP_dd = r1 * (cos(phi_i) * phidot_i**2 + sin(phi_i) * ddphi)
    yP_dd = r1 * (sin(phi_i) * phidot_i**2 - cos(phi_i) * ddphi)
    xB_dd = xP_dd - l * sin(theta_i) * thetadot_i**2 + l * cos(theta_i) * ddtheta
    yB_dd = yP_dd + l * cos(theta_i) * thetadot_i**2 + l * sin(theta_i) * ddtheta
    
    e_r = np.array([xP[i] - xB[i], yP[i] - yB[i]])
    e_r /= np.linalg.norm(e_r)
    T_bob[i] = e_r.dot(m * np.array([xB_dd, yB_dd]) + np.array([0, m*g]))

    xA_dd = -r2 * (cos(phi_i) * phidot_i**2 + sin(phi_i) * ddphi)
    yA_dd = -r2 * (sin(phi_i) * phidot_i**2 - cos(phi_i) * ddphi)
    T_M[i] = M * (yA_dd + g)
    
    F_rod_bob = -T_bob[i] * e_r
    rP = np.array([xP[i], yP[i]])
    tau_bob[i] = rP[0]*F_rod_bob[1] - rP[1]*F_rod_bob[0]
    
    F_rod_M = np.array([0, -T_M[i]])
    rA = np.array([xA[i], yA[i]])
    tau_M[i] = rA[0]*F_rod_M[1] - rA[1]*F_rod_M[0]

# ============================
# VISUALIZACIÓN SINCRONIZADA
# ============================
fig = plt.figure(figsize=(10,10))
gs = fig.add_gridspec(4, 1, height_ratios=[2,1,1,1])
ax_anim = fig.add_subplot(gs[0])
ax_ang  = fig.add_subplot(gs[1])
ax_forc = fig.add_subplot(gs[2])
ax_en   = fig.add_subplot(gs[3])
plt.subplots_adjust(hspace=0.45)

# --- Animación ---
ax_anim.set_aspect('equal')
Lmax = r1 + r2 + l + 0.5
ax_anim.set_xlim(-Lmax, Lmax)
ax_anim.set_ylim(-l - 1.0, Lmax * 0.2)
ax_anim.set_title('Sistema acoplado sin amortiguamiento (energía conservada)')
ax_anim.axis('off')

rod_line, = ax_anim.plot([], [], 'k-', lw=4)
pend_line, = ax_anim.plot([], [], color='tab:blue', lw=2)
bob_dot, = ax_anim.plot([], [], 'o', color='tab:blue', ms=8)
mass_line, = ax_anim.plot([], [], color='k', lw=2)
mass_dot, = ax_anim.plot([], [], 's', color='dimgray', ms=10)
fulcrum, = ax_anim.plot(0, 0, 'k^', ms=10)

# --- Subgráficos ---
# Ángulos
ax_ang.set_xlim(0, t_max)
ax_ang.set_ylim(-0.2, 0.2)
ax_ang.set_xlabel('t (s)')
ax_ang.set_ylabel('Ángulo [rad]')
ax_ang.set_title('Ángulos')
line_phi, = ax_ang.plot([], [], label=r'$\varphi$ (varilla)')
line_theta, = ax_ang.plot([], [], label=r'$\theta$ (péndulo)')
ax_ang.legend(); ax_ang.grid(True)

# Tensiones y torques
ax_forc.set_xlim(0, t_max)
ax_forc.set_ylim(-10, 10)
ax_forc.set_xlabel('t (s)')
ax_forc.set_ylabel('T [N],  τ [N·m]')
ax_forc.set_title('Tensiones y torques')
line_Tb, = ax_forc.plot([], [], color='tab:blue', label=r'$T_{bob}$ [N]')
line_Tm, = ax_forc.plot([], [], color='tab:purple', label=r'$T_M$ [N]')
line_tb, = ax_forc.plot([], [], color='tab:orange', label=r'$\tau_{bob}$ [N·m]')
line_tm, = ax_forc.plot([], [], color='tab:green', label=r'$\tau_M$ [N·m]')
ax_forc.legend(); ax_forc.grid(True)

# Energía total
ax_en.set_xlim(0, t_max)
E0 = E_total[0]
ax_en.set_ylim(E0 - 2, E0 + 2)
ax_en.set_xlabel('t (s)')
ax_en.set_ylabel('Energía [J]')
ax_en.set_title('Energía total del sistema')
line_E, = ax_en.plot([], [], color='tab:red', label='Energía total [J]')
ax_en.legend(); ax_en.grid(True)

# ============================
# ANIMACIÓN SINCRONIZADA
# ============================
skip = 20
frames = range(0, len(t), skip)

def init():
    for line in [rod_line, pend_line, bob_dot, mass_line, mass_dot,
                 line_phi, line_theta, line_Tb, line_Tm, line_tb, line_tm, line_E]:
        line.set_data([], [])
    return [rod_line, pend_line, bob_dot, mass_line, mass_dot,
            line_phi, line_theta, line_Tb, line_Tm, line_tb, line_tm, line_E]

def update(i):
    # Movimiento
    rod_line.set_data([xP[i], xA[i]], [yP[i], yA[i]])
    pend_line.set_data([xP[i], xB[i]], [yP[i], yB[i]])
    bob_dot.set_data(xB[i], yB[i])
    mass_line.set_data([xA[i], xM[i]], [yA[i], yM[i]])
    mass_dot.set_data(xM[i], yM[i])

    # Gráficos
    tt = t[:i]
    line_phi.set_data(tt, phi[:i])
    line_theta.set_data(tt, theta[:i])
    line_Tb.set_data(tt, T_bob[:i])
    line_Tm.set_data(tt, T_M[:i])
    line_tb.set_data(tt, tau_bob[:i])
    line_tm.set_data(tt, tau_M[:i])
    line_E.set_data(tt, E_total[:i])
    return [rod_line, pend_line, bob_dot, mass_line, mass_dot,
            line_phi, line_theta, line_Tb, line_Tm, line_tb, line_tm, line_E]

ani = FuncAnimation(fig, update, frames=frames, init_func=init, blit=True, interval=30)
plt.show()
