import numpy as np

# Parameters
L = 200e-6
epsilon = 0.8
alpha = 2
sigma = 1e-5
dt = 0.1
m_cell = 0.1
tol = 1e-8
lim = [0, L]
m = 9
n = 3

# Sampling
x0 = 1e-8  # Start from 1e-8 to avoid division by zero
xf = L
n_samples = 100
dx = (xf - x0) / n_samples

# Compute LJ force
x_values = []
lj_forces = []
xi = x0

while xi <= xf:
    force = alpha * epsilon * (
        m * sigma**m / xi**(m + 1)
        - n * sigma**n / xi**(n + 1)
    )
    x_values.append(xi)
    lj_forces.append(force)
    xi += dx

# Save to TXT
# Save to TXT (single line, whitespace-separated)
with open("lj_force_data.txt", "w") as f:
    f.write(" ".join(f"{f_val:.8e}" for f_val in lj_forces))
