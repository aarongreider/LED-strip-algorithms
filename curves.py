import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import random

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
    """ for p0, p1, p2 in curves:
        curve_pts = [quad_bezier_point(p0, p1, p2, t/1000) for t in range(1001)]
        cx, cy = zip(*curve_pts)
        plt.plot(cx, cy, 'b-', alpha=0.3) """

    # Plot the sampled points spaced approximately equally
    px, py = zip(*points)
    plt.plot(px, py, 'o', label='Evenly spaced points')

    # Mark control points
    """ for seg in curves:
            for pt in seg:
                plt.plot(pt[0], pt[1], 'go') """

    plt.title("Stitched Quadratic Bézier Curves with Equal Spacing")
    plt.axis('equal')
    plt.legend()
    plt.show()


if __name__ == "__main__":
    # Define multiple quadratic Bézier curves forming a loop
    curves = [
        ((0, 0), (10, 0), (10, 10)),
        ((10, 10), (6, 14), (2, 10)),
        ((2, 10), (2, 0), (12, 0)),
    ]

    steps_per_curve = 200

    # Get evenly spaced points along all curves
    points = stitch_quadratic_beziers(curves, 144, steps_per_curve)
    plot_points(curves, points)
