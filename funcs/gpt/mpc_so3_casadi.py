# ============================================
# MPC for Attitude Control on SO(3) using CasADi
# ============================================
import casadi as ca
import numpy as np
from scipy.linalg import expm

# Parameters
J = np.diag([1.0, 1.2, 1.5])
Jinv = np.linalg.inv(J)
dt = 0.05
N = 20  # Horizon

# CasADi symbols
omega = ca.SX.sym('omega', 3)
R = ca.SX.sym('R', 3, 3)
tau = ca.SX.sym('tau', 3)
x = ca.vertcat(ca.reshape(R, 9, 1), omega)


def rodrigues(omega):
    theta = ca.norm_2(omega)
    omega_hat = hat(omega)
    I = ca.MX.eye(3)

    small = 1e-6

    A = ca.if_else(theta > small, ca.sin(theta)/theta, 1 - theta**2 / 6)
    B = ca.if_else(theta > small, (1 - ca.cos(theta)) / (theta**2), 0.5 - theta**2 / 24)

    return I + A * omega_hat + B * (omega_hat @ omega_hat)



def hat(w):
    return ca.vertcat(
        ca.horzcat(0, -w[2], w[1]),
        ca.horzcat(w[2], 0, -w[0]),
        ca.horzcat(-w[1], w[0], 0)
    )


# Dynamics: x_next = f(x, u)
def dynamics(x, u):
    R = ca.reshape(x[0:9], 3, 3)
    omega = x[9:12]
    # R_next = R @ ca.expm(dt * hat(omega))
    # R_next = Rk @ expm(dt * hat(omegak))
    R_next = R @ rodrigues(dt * omega)

    omega_dot = Jinv @ (u - ca.cross(omega, J @ omega))
    omega_next = omega + dt * omega_dot
    return ca.vertcat(ca.reshape(R_next, 9, 1), omega_next)


# Build MPC problem
opti = ca.Opti()
X = opti.variable(12, N + 1)  # State trajectory
U = opti.variable(3, N)  # Control inputs

x0_param = opti.parameter(12)  # Initial state
Rd_param = opti.parameter(3, 3)  # Desired rotation matrix

# Objective function
cost = 0
for k in range(N):
    Rk = ca.reshape(X[0:9, k], 3, 3)
    omegak = X[9:12, k]
    Re = Rd_param.T @ Rk
    phi = 0.5 * (Re - Re.T)
    err_vec = ca.vertcat(phi[2, 1], phi[0, 2], phi[1, 0])  # vee(log approx)
    cost += ca.mtimes([err_vec.T, err_vec]) + 0.1 * ca.mtimes([U[:, k].T, U[:, k]])

    # Dynamics constraints
    x_next = dynamics(X[:, k], U[:, k])
    opti.subject_to(X[:, k + 1] == x_next)

# Initial condition constraint
opti.subject_to(X[:, 0] == x0_param)

# Input bounds
opti.subject_to(opti.bounded(-5, U, 5))

# Solver
opti.minimize(cost)
opti.solver('ipopt')

# Simulation loop
x_sim = np.zeros((12, 100))
omega0 = np.array([0.5, -0.2, 0.3])
R0 = expm(hat(np.array([0.3, 0.2, 0.1])))
x0 = np.concatenate((R0.flatten(), omega0))
x_sim[:, 0] = x0
Rd = np.eye(3)

for t in range(99):
    opti.set_value(x0_param, x0)
    opti.set_value(Rd_param, Rd)

    sol = opti.solve()
    u_opt = sol.value(U[:, 0])

    # Integrate dynamics forward
    Rk = np.reshape(x0[:9], (3, 3))
    omegak = x0[9:12]
    R_next = Rk @ expm(dt * hat(omegak))
    omega_dot = Jinv @ (u_opt - np.cross(omegak, J @ omegak))
    omega_next = omegak + dt * omega_dot
    x0 = np.concatenate((R_next.flatten(), omega_next))

    x_sim[:, t + 1] = x0

# Plotting (example: omega)
import matplotlib.pyplot as plt

plt.plot(np.linspace(0, dt * 100, 100), x_sim[9:12, :].T)
plt.xlabel('Time [s]')
plt.ylabel('Angular Velocity [rad/s]')
plt.title('Angular Velocity under MPC')
plt.legend(['ω₁', 'ω₂', 'ω₃'])
plt.grid()
plt.show()
