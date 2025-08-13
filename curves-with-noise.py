import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import random
import math


def quad_bezier_point(p0, p1, p2, t):
    """Quadratic Bézier point at t (0..1)."""
    u = 1 - t
    return (
        u*u*p0[0] + 2*u*t*p1[0] + t*t*p2[0],
        u*u*p0[1] + 2*u*t*p1[1] + t*t*p2[1]
    )


def approx_bezier_length(p0, p1, p2, steps=20):
    """Approximate length of a quadratic Bézier by sampling."""
    length = 0
    prev = quad_bezier_point(p0, p1, p2, 0)
    for i in range(1, steps+1):
        t = i / steps
        curr = quad_bezier_point(p0, p1, p2, t)
        dx = curr[0] - prev[0]
        dy = curr[1] - prev[1]
        length += (dx*dx + dy*dy) ** 0.5
        prev = curr
    return length


def stitch_quadratic_beziers(curves, num_points, steps_per_curve=100):
    """
    Generate exactly num_points evenly spaced along stitched quadratic Bézier curves.

    curves: list of (p0, p1, p2) control point triples.
    num_points: total number of points desired along entire stitched path.
    steps_per_curve: samples per segment for length estimation and stepping.

    Returns: list of (x, y) points.
    """
    # Precompute lengths of each curve
    lengths = []
    total_length = 0
    for p0, p1, p2 in curves:
        seg_len = approx_bezier_length(p0, p1, p2, steps=steps_per_curve)
        lengths.append(seg_len)
        total_length += seg_len

    spacing = total_length / (num_points - 1)  # distance between points

    points = []
    seg_index = 0
    dist_in_seg = 0

    p0, p1, p2 = curves[0]
    prev_point = quad_bezier_point(p0, p1, p2, 0)
    points.append(prev_point)

    dist_covered = 0  # total distance walked so far

    # To keep track of the t parameter on the current segment
    t_index = 0

    while len(points) < num_points and seg_index < len(curves):
        p0, p1, p2 = curves[seg_index]

        while t_index < steps_per_curve:
            t_index += 1
            t = t_index / steps_per_curve
            curr_point = quad_bezier_point(p0, p1, p2, t)
            dx = curr_point[0] - prev_point[0]
            dy = curr_point[1] - prev_point[1]
            step_len = (dx*dx + dy*dy) ** 0.5
            dist_in_seg += step_len
            dist_covered += step_len

            if dist_in_seg >= spacing:
                points.append(curr_point)
                dist_in_seg = 0
                prev_point = curr_point
                if len(points) >= num_points:
                    break
            else:
                prev_point = curr_point

        # Move to next segment
        seg_index += 1
        t_index = 0
        dist_in_seg = 0
        if seg_index < len(curves):
            prev_point = quad_bezier_point(*curves[seg_index], 0)

    # If somehow points are less than requested (due to rounding), pad with last point
    while len(points) < num_points:
        points.append(points[-1])

    return points


def plot_points(curves, points):
    """ Plot the full smooth curve for each segment (fine sampling) """
    plt.figure(figsize=(8, 8))

    # Plot the sampled points spaced approximately equally
    px, py = zip(*points)

    plt.plot(px, py, 'o', label='Evenly spaced points')
    plt.title("Stitched Quadratic Bézier Curves with Equal Spacing")
    plt.axis('equal')
    plt.legend()
    plt.show()


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
    n = (a * (1 - u[0]) + b * u[0]) * (1 - u[1]) + \
           (c * (1 - u[0]) + d * u[0]) * u[1]
    
    # Map noise value [0,1) to palette index
    """ idx = int(n * len(color_palette)) % len(color_palette)
    
    return color_palette[idx] """
    return n

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

# Define multiple quadratic Bézier curves forming a loop
curves = [
    ((0, 0), (10, 0), (10, 10)),
    ((10, 10), (6, 14), (2, 10)),
    ((2, 10), (2, 0), (12, 0)),
]
steps_per_curve = 200  # sampling steps
num_points = 144

# Get evenly spaced points along all curves
points = stitch_quadratic_beziers(curves, num_points, steps_per_curve)

# Plotting
fig, ax = plt.subplots()
px, py = zip(*points)  # px = all x values, py = all y values
dots = ax.scatter(px, py, c='r', s=80)
ani = FuncAnimation(fig, update, frames=200, interval=50, blit=True)
plt.show()

# plot_points(curves, points)
