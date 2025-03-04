import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Define the polynomial p(x) = x^3 - 1 and its evaluation
def poly_eval(z, coeffs):
    return np.polyval(coeffs, z)

# Durand-Kerner one iteration update
def durand_kerner_step(z, coeffs):
    new_z = np.copy(z)
    for i in range(len(z)):
        prod = 1
        for j in range(len(z)):
            if i != j:
                prod *= (z[i] - z[j])
        new_z[i] = z[i] - poly_eval(z[i], coeffs) / prod
    return new_z

# Parameters: polynomial x^3 - 1
coeffs = [1, 0, 0, -1]  # Represents x^3 - 1
true_roots = np.roots(coeffs)

# Initialize the guesses (shifted away from the true roots)
n_roots = 3
radius = 1.2  # Different from the unit circle
angles = np.linspace(np.pi, 2*np.pi, n_roots, endpoint=False)
z = radius * np.exp(1j * angles)

# Precompute iterations for the animation
iterations = [z.copy()]
max_iters = 30
tol = 1e-8

for _ in range(max_iters):
    z_next = durand_kerner_step(z, coeffs)
    iterations.append(z_next.copy())
    if np.all(np.abs(z_next - z) < tol):
        break
    z = z_next

iterations = np.array(iterations)

# Set up the plot for animation
fig, ax = plt.subplots(figsize=(6,6))
ax.set_title("Durand-Kerner Algorithm Convergence")
ax.set_xlabel("Real part")
ax.set_ylabel("Imaginary part")
ax.grid(True)
ax.set_aspect('equal')

# Plot the true roots for reference
ax.scatter(true_roots.real, true_roots.imag, color='red', marker='x', s=100, label='True Roots')
ax.legend()

ax.set_xlim(-2, 2)
ax.set_ylim(-2, 2)

scat = ax.scatter([], [], color='blue', s=100, label='Current Estimates')

def init():
    scat.set_offsets(np.empty((0, 2)))
    return scat,

def update(frame):
    current_estimates = iterations[frame]
    points = np.column_stack((current_estimates.real, current_estimates.imag))
    scat.set_offsets(points)
    ax.set_title(f"Iteration {frame}")
    return scat,

ani = animation.FuncAnimation(fig, update, frames=len(iterations),
                              init_func=init, blit=True, interval=500, repeat=False)

plt.show()
