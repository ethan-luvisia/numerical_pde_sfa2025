# %%
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import solve_banded
from mpl_toolkits.mplot3d import Axes3D

# testing factor
S_max = 150
K = 100
T = 1.0
r = 0.05
sigma0 = 0.2
alpha = 0.1
beta = 0.05

M = 300          # space grid
N = 2000         # time grid

dS = S_max / M
dt = T / N

S = np.linspace(0, S_max, M + 1)
t = np.linspace(0, T, N + 1)

# initial condition
V = np.maximum(S - K, 0)
solution = np.zeros((N + 1, M + 1))
solution[-1, :] = V.copy()

# boundary condition
def boundary_conditions(tau):
    left = 0.0
    right = S_max - K * np.exp(-r * (T - tau))
    return left, right

# local vol func
def local_vol(S, tau):
    sig = sigma0 * (1 + alpha * np.sin(np.pi * S / K) + beta * tau)
    return np.clip(sig, 0.05, 0.4)

# coeffs
def compute_coeffs(Si, sig, dS, r, dt):
    Si2 = np.clip(Si ** 2, 1e-8, 1e6)
    A = 0.25 * dt * (sig ** 2 * Si2 / dS ** 2 - r * Si / dS)
    B = -0.5 * dt * (sig ** 2 * Si2 / dS ** 2 + r)
    C = 0.25 * dt * (sig ** 2 * Si2 / dS ** 2 + r * Si / dS)
    return A, B, C

# Crank-Nicolson main loop
for n in reversed(range(N)):
    tau = t[n]
    sigma_vec = local_vol(S, tau)

    if np.any(np.isnan(solution[n+1])) or np.any(np.isinf(solution[n+1])):
        raise ValueError(f"Numerical instability at step {n+1}")

    a = np.zeros(M - 1)
    b = np.zeros(M - 1)
    c = np.zeros(M - 1)
    rhs = np.zeros(M - 1)

    # factor matrix
    for i in range(1, M):
        Si = S[i]
        sig = sigma_vec[i]
        A, B, C = compute_coeffs(Si, sig, dS, r, dt)
        a[i - 1] = -A
        b[i - 1] = 1 - B
        c[i - 1] = -C

    # right vector
    for i in range(1, M):
        Si = S[i]
        sig = sigma_vec[i]
        A, B, C = compute_coeffs(Si, sig, dS, r, -dt)  # backward
        rhs[i - 1] = A * solution[n + 1, i - 1] + (1 + B) * solution[n + 1, i] + C * solution[n + 1, i + 1]

    # boundary condition
    left, right = boundary_conditions(tau)
    rhs[0] -= a[0] * left
    rhs[-1] -= c[-1] * right

    # solving linear system
    ab = np.zeros((3, M - 1))
    ab[0, 1:] = c[:-1]
    ab[1, :] = b
    ab[2, :-1] = a[1:]
    V_inner = solve_banded((1, 1), ab, rhs)

    solution[n, 0] = left
    solution[n, M] = right
    solution[n, 1:M] = V_inner

# graph
T_grid, S_grid = np.meshgrid(t, S)

fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(S_grid, T_grid, solution.T, cmap='viridis')
ax.set_xlabel('Stock Price S')
ax.set_ylabel('Time t')
ax.set_zlabel('Option Price V')
ax.set_title('European Call Option under Local Volatility (Stable)')
plt.tight_layout()
plt.show()


# %%
from scipy.interpolate import UnivariateSpline
from mpl_toolkits.mplot3d import Axes3D

# initial spline
delta_spline = np.zeros_like(solution)
gamma_spline = np.zeros_like(solution)

# for each time layer t[n]，compute spline and take derivative
for n in range(N + 1):
    spline = UnivariateSpline(S, solution[n], s=0.5)  # s=0.5 control sommth
    delta_spline[n] = spline.derivative(n=1)(S)
    gamma_spline[n] = spline.derivative(n=2)(S)

# clip get rid of max
delta_spline = np.clip(delta_spline, -2, 2)
gamma_spline = np.clip(gamma_spline, -0.2, 0.5)

# create grid
T_grid, S_grid = np.meshgrid(t, S)

# draw Delta Surface
fig1 = plt.figure(figsize=(10, 7))
ax1 = fig1.add_subplot(111, projection='3d')
surf1 = ax1.plot_surface(S_grid, T_grid, delta_spline.T, cmap='viridis')
ax1.set_xlabel('Stock Price S')
ax1.set_ylabel('Time t')
ax1.set_zlabel('Delta')
ax1.set_title('Smoothed Delta Surface (Spline Fit + Softplus)')
fig1.colorbar(surf1, ax=ax1, shrink=0.6, aspect=10)
plt.tight_layout()
plt.show()
