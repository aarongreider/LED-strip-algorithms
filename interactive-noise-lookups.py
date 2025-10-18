import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import math
import time
import colorsys
import random




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
    """Linear interpolation between a and b by normalized value t."""
    return a * (1 - t) + b * t


def inverse_lerp(a, b, x):
    """ Inverse Interpolation to normalize x by a and b"""
    return (x - a) / (b - a)


def nudge(current, target, up_speed, down_speed):
    """ Eases or nudges the current toward the target by factor speed """
    if current < target:
        return min(current + up_speed, target)
    elif current > target:
        return max(current - down_speed, target)
    return current


def height_to_hue(height, hues):
    """
    Map a normalized value [0,1] to a color from a list of hue integers.
    Smoothly interpolates between neighboring hues.
    """
    height = clamp(height)
    i, j, t = get_color_indices(height, len(hues))
    j = min(j, len(hues) - 1)  # ensure j doesnt go out of bounds
    return lerp(hues[i],  hues[j], t)

def get_height(x, y, scale):
    """ Gets the intensity of the noise function based on position x, y, and scale """
    """
    x, y: floats, scaled position
    Returns: interpolated noise value 0-1
    """
    x = x * scale
    y = y * scale
    # Wrap around the grid
    gx = int(x) % GRID_SIZE
    gy = int(y) % GRID_SIZE
    gx1 = (gx + 1) % GRID_SIZE
    gy1 = (gy + 1) % GRID_SIZE

    fx = x - int(x)  # fractional part
    fy = y - int(y)

    # Four corners
    a = noise_grid[gy][gx]
    b = noise_grid[gy][gx1]
    c = noise_grid[gy1][gx]
    d = noise_grid[gy1][gx1]

    # Linear interpolation
    top = a + (b - a) * fx
    bottom = c + (d - c) * fx
    return top + (bottom - top) * fy


def get_distance():
    """ generates a sin betwixt 0 and 10 if distance, the unit is meters """
    # return (math.sin(t/2) + 1) * 5
    return mouse_x


def refresh_values():
    # nudge values toward target
    # set new targets
    for i, value in enumerate(values):
        values[i] = nudge(value, target_values[i], flutter_speed_up, flutter_speed_down)
        if random.uniform(0, 1) < flutter_probability: # add random flutters
            target_values[i] = 1.0
        if value == target_values[i]:
            target_values[i] = min_v

offset = 0
def set_hsv(frame):
    """ Set the LEDS """
    global flutter_probability, flutter_speed_up, flutter_speed_down, min_v
    global offset

    
    scale = clamp(inverse_lerp(-11, 11, get_distance()), 0, 1) # normalized distance

    # smaller scale → faster x, y offset
    min_speed = 0.002
    max_speed = .5
    speed = min_speed + (1 - scale) * (max_speed - min_speed)

    # Increment offset linearly, scaled by speed
    offset += speed
    print('%.3f'%speed, '%.3f'%offset)

    flutter_probability = lerp(.2, .02, scale)
    flutter_speed_up = lerp(.15, .1, scale)
    flutter_speed_down = lerp(.09, .02, scale)
    min_v = lerp(.75, .25, scale)
   #  print("flutter:", flutter_probability, flutter_speed_up, flutter_speed_down)

    refresh_values()

    for (x, y), i in points:
        # Add time offset for flowing noise3
        # print("getting height with", x, y, offset, scale)
        height = get_height(x + offset, y + offset * 0.3, scale)
        # Map noise to leds
        hue = height_to_hue(height, hue_palette)
        # normalize and invert value for hsv
        saturation = lerp(.75, 1, 1 - scale)
        hsv = (hue, saturation, values[i])
        # print("scales", scale, sv_scale)
        # print('hsv', hsv)
        LEDS[i] = (hsv)

    # set dots with converted hsv
    dots.set_facecolors([colorsys.hsv_to_rgb(*hsv) for hsv in LEDS])
    time.sleep(.01)
    return [dots]


NUM_LEDS = 120
LEDS = [(0, 0, 0)] * NUM_LEDS  # preallocate memory
print(LEDS)

flutter_probability = .02
flutter_speed_up = .05
flutter_speed_down = .2
min_v = .38
values = [min_v] * NUM_LEDS
target_values = [min_v] * NUM_LEDS # the target to lerp towards


# Get evenly spaced points along all curves
points = get_points_on_circle((0, 0), 10, NUM_LEDS)
print(points)

hue_palette = [0, .05, .1, .2, .3]
# hue_palette = [.3, .38, .4, .58, .62]
GRID_SIZE = 16
noise_grid = [[random.random() for _ in range(GRID_SIZE)] for _ in range(GRID_SIZE)]


# region replace plotting callbacks with LED loop
# Plotting
fig, ax = plt.subplots()
px, py = zip(*[p for p, _ in points])  # px = all x values, py = all y values
dots = ax.scatter(px, py, c='r', s=80)

mouse_x = 0
mouse_y = 0

def on_mouse_move(event):
    global mouse_x, mouse_y
    if event.inaxes:  # Only track when inside the plot
        mouse_x, mouse_y = event.xdata, event.ydata
        # print("Mouse at data coords:", mouse_x, mouse_y)


# Connect the motion event to the handler
cid = fig.canvas.mpl_connect('motion_notify_event', on_mouse_move)


ani = FuncAnimation(fig, set_hsv, frames=1000, interval=50, blit=False)
plt.show()
# endregion