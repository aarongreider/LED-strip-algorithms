import math
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# --- Random function ---
def random2(st):
    """2D random function like GLSL."""
    dot = st[0] * 12.9898 + st[1] * 78.233
    return (math.sin(dot) * 43758.5453123) % 1.0

# --- 2D noise function ---
def noise2(st):
    """2D noise using bilinear interpolation."""
    # Floor and fract
    i = [math.floor(st[0]), math.floor(st[1])]
    f = [st[0] - i[0], st[1] - i[1]]

    # Four corners
    a = random2(i)
    b = random2([i[0] + 1.0, i[1]])
    c = random2([i[0], i[1] + 1.0])
    d = random2([i[0] + 1.0, i[1] + 1.0])

    # Smoothstep (cubic Hermite)
    u = [f[0]*f[0]*(3.0 - 2.0*f[0]), f[1]*f[1]*(3.0 - 2.0*f[1])]

    # Bilinear interpolation
    return (a * (1 - u[0]) + b * u[0]) * (1 - u[1]) + \
           (c * (1 - u[0]) + d * u[0]) * u[1]

# --- Data points ---
num_points = 144
px = [2 * math.pi * i / (num_points - 1) for i in range(num_points)]
py = [math.sin(x) for x in px]

fig, ax = plt.subplots()
dots = ax.scatter(px, py, c='r', s=80)

# Pack points as list of [x, y]
points = [[x, y] for x, y in zip(px, py)]

def update(frame):
    t = frame * 0.05
    colors = []
    for p in points:
        # Add time offset for flowing noise
        n = noise2([p[0] + t, p[1] + t * 0.3])
        # Map noise to RGB
        colors.append([n, 1 - n, 0.5])
    dots.set_facecolors(colors)
    return dots,

ani = FuncAnimation(fig, update, frames=200, interval=50, blit=True)
plt.show()
