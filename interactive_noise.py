import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import math
import time


def get_points_on_circle(center, radius, num_points):
    cx, cy = center
    return [
        (cx + radius * math.cos(2*math.pi*i/num_points),
         cy + radius * math.sin(2*math.pi*i/num_points))
        for i in range(num_points)
    ]

def smooth_random(ix, iy):
    # deterministic non-random function, accounts for position of value
    val = math.sin(ix * 12.9898 + iy * 78.233) * 43758.5453 % 1
    print("random value", val, "at", ix, iy)
    return val

def noise2(x, y, scale):
    x = x * scale
    y = y * scale
    # Bilinear Interpolation function   
    cell_x = math.floor(x)   # -> X coordinate (lower-left corner)
    cell_y = math.floor(y)   # -> Y coordinate (lower-left corner)
    offset_x = x - cell_x    # -> fractional offset inside the cell along X
    offset_y = y - cell_y    # -> fractional offset inside the cell along Y
    # pass in the boundaries
    a = smooth_random(cell_x, cell_y)
    b = smooth_random(cell_x + 1.0, cell_y)
    c = smooth_random(cell_x, cell_y + 1.0)
    d = smooth_random(cell_x + 1.0, cell_y + 1.0)   
    u = offset_x * offset_x * (3.0 - 2.0 * offset_x)
    v = offset_y * offset_y * (3.0 - 2.0 * offset_y)
    return (a * (1 - u) + b * u) * (1 - v) + (c * (1 - u) + d * u) * v

def get_distance(t):
    # generates a sin betwixt 0 and 10
    # if distance, the unit is meters
    return (math.sin(t/2) + 1) * 5

def update(frame):
    time = frame * 0.05
    scale = get_distance(time) / 5
    colors = []
    for point in points:
        # Add time offset for flowing noise3
        noise = noise2(point[0] + time, point[1] + time * 0.3, scale)
        # Map noise to RGB #TODO: make hsv
        colors.append([noise, 1 - noise, 0.5])
    dots.set_facecolors(colors)
    return [dots]


NUM_LEDS = 120

# Get evenly spaced points along all curves
points = get_points_on_circle((0,0), 10, NUM_LEDS)

# Plotting
fig, ax = plt.subplots()
px, py = zip(*points)  # px = all x values, py = all y values
dots = ax.scatter(px, py, c='r', s=80)
ani = FuncAnimation(fig, update, frames=200, interval=50, blit=True)
plt.show()
time.sleep(.2)