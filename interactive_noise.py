import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import math
import time
import colorsys

# region


def get_points_on_circle(center, radius, num_points):
    cx, cy = center
    return [
        ((
            cx + radius * math.cos(2 * math.pi * i / num_points),
            cy + radius * math.sin(2 * math.pi * i / num_points)
        ), i)
        for i in range(num_points)
    ]


def clamp(value, min_val=0.0, max_val=1.0):
    """Clamp a number between min_val and max_val."""
    return max(min_val, min(max_val, value))


def get_color_indices(value, num_colors):
    """
    Given a normalized value [0,1], find the two color indices
    to interpolate between and the local interpolation factor t.
    """
    n = num_colors - 1
    scaled = value * n
    i = int(scaled)
    j = i + 1
    t = scaled - i
    return i, j, t


def lerp(a, b, t):
    """Linear interpolation between a and b by t."""
    return a * (1 - t) + b * t


def height_to_hue(height, hue_palette):
    """
    Map a normalized value [0,1] to a color from a list of hue integers.
    Smoothly interpolates between neighboring hues.
    """
    height = clamp(height)
    i, j, t = get_color_indices(height, len(hue_palette))
    j = min(j, len(hue_palette) - 1)  # ensure j doesnt go out of bounds
    return lerp(hue_palette[i],  hue_palette[j], t)


def get_random_value(ix, iy):
    # deterministic random function, accounts for position of value
    val = math.sin(ix * 12.9898 + iy * 78.233) * 43758.5453
    val = val - math.floor(val)
    return val


def get_height(x, y, scale):
    """ Gets the intensity of the noise function based on position x, y, and scale """
    x = x * scale
    y = y * scale
    # Bilinear Interpolation function
    cell_x = math.floor(x)   # -> X coordinate (lower-left corner)
    cell_y = math.floor(y)   # -> Y coordinate (lower-left corner)
    offset_x = x - cell_x    # -> fractional offset inside the cell along X
    offset_y = y - cell_y    # -> fractional offset inside the cell along Y
    # pass in the boundaries
    a = get_random_value(cell_x, cell_y)
    b = get_random_value(cell_x + 1.0, cell_y)
    c = get_random_value(cell_x, cell_y + 1.0)
    d = get_random_value(cell_x + 1.0, cell_y + 1.0)
    u = offset_x * offset_x * (3.0 - 2.0 * offset_x)  # smoothstep function
    v = offset_y * offset_y * (3.0 - 2.0 * offset_y)  # smoothstep function

    instensity = (a * (1 - u) + b * u) * (1 - v) + (c * (1 - u) + d * u) * v
    return instensity


def get_distance(t):
    # generates a sin betwixt 0 and 10
    # if distance, the unit is meters
    return (math.sin(t/2) + 1) * 5


def set_hsv(frame):
    """ Set the LEDS """
    offset = frame * 0.05
    scale = get_distance(offset) / 5

    for (x, y), i in points:
        # Add time offset for flowing noise3
        print("getting height with", x, y, offset, scale)
        height = get_height(x + offset, y + offset * 0.3, scale)
        # Map noise to leds
        hue = height_to_hue(height, hue_palette)
        hsv = (hue, 1.0, 1.0)
        print('hsv', hsv)
        LEDS[i] = (hsv)

    # set dots with converted hsv
    dots.set_facecolors([colorsys.hsv_to_rgb(*hsv) for hsv in LEDS])
    return [dots]

print("Hello World")
NUM_LEDS = 120
LEDS = [(0, 0, 0)] * NUM_LEDS  # preallocate memory
print(LEDS)

# Get evenly spaced points along all curves
points = get_points_on_circle((0, 0), 10, NUM_LEDS)
print(points)

hue_palette = [0, .25, .5, .75, 1.0]

# Plotting
fig, ax = plt.subplots()
px, py = zip(*[p for p, _ in points])  # px = all x values, py = all y values
dots = ax.scatter(px, py, c='r', s=80)
ani = FuncAnimation(fig, set_hsv, frames=200, interval=50, blit=False)
plt.show()
time.sleep(.2)
