import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# --- Noise functions from GLSL ---
def random2(st):
    dot = st[..., 0] * 12.9898 + st[..., 1] * 78.233
    return np.mod(np.sin(dot) * 43758.5453123, 1.0)

def noise2(st):
    i = np.floor(st)
    f = st - i  # fract
    a = random2(i)
    b = random2(i + [1.0, 0.0])
    c = random2(i + [0.0, 1.0])
    d = random2(i + [1.0, 1.0])
    u = f * f * (3.0 - 2.0 * f)
    return (a * (1 - u[..., 0]) + b * u[..., 0]) * (1 - u[..., 1]) + \
           (c * (1 - u[..., 0]) + d * u[..., 0]) * u[..., 1]

# --- Data points ---
px = np.linspace(0, 2*np.pi, 10)
py = np.sin(px)

fig, ax = plt.subplots()
dots = ax.scatter(px, py, c='r', s=80)

# Normalize coords for noise lookup
points = np.stack([px, py], axis=-1)

def update(frame):
    # Add "time" to the noise input
    t = frame * 0.05
    n_vals = noise2(points + [t, t])
    # Map noise (0-1) to RGB smoothly
    colors = np.stack([n_vals, 1 - n_vals, 0.5 * np.ones_like(n_vals)], axis=-1)
    dots.set_facecolors(colors)
    return dots,

ani = FuncAnimation(fig, update, frames=200, interval=50, blit=True)
plt.show()
