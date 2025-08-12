import matplotlib.pyplot as plt

def quad_bezier_point(p0, p1, p2, t):
    """Quadratic Bézier equation."""
    u = 1 - t
    return (
        u*u*p0[0] + 2*u*t*p1[0] + t*t*p2[0],
        u*u*p0[1] + 2*u*t*p1[1] + t*t*p2[1]
    )

def bezier_equal_spacing_gen(p0, p1, p2, spacing, steps=1000):
    """Yield approximately equally spaced points along a quadratic Bézier."""
    prev = quad_bezier_point(p0, p1, p2, 0)
    yield prev
    dist_acc = 0
    for i in range(1, steps + 1):
        t = i / steps
        curr = quad_bezier_point(p0, p1, p2, t)

        dx = curr[0] - prev[0]
        dy = curr[1] - prev[1]
        seg_len = (dx*dx + dy*dy) ** 0.5  # Euclidean distance

        dist_acc += seg_len
        if dist_acc >= spacing:
            yield curr
            dist_acc = 0
        prev = curr


if __name__ == "__main__":
    p0 = (0, 0)
    p1 = (5, 10)
    p2 = (10, 0)
    spacing = 1.0

    # Generate points with equal arc length
    spaced_points = list(bezier_equal_spacing_gen(p0, p1, p2, spacing=spacing, steps=500))

    # Generate a smooth curve for reference
    curve_t = [quad_bezier_point(p0, p1, p2, t/1000) for t in range(1001)]
    curve_x, curve_y = zip(*curve_t)

    # Extract spaced points into x/y lists
    px, py = zip(*spaced_points)

    # Plot
    plt.figure(figsize=(6, 6))
    plt.plot(curve_x, curve_y, 'b-', label='Bézier curve')
    plt.plot(px, py, 'ro', label=f'Points ~{spacing} units apart')
    plt.plot(*p0, 'go', label='P0')
    plt.plot(*p1, 'yo', label='P1')
    plt.plot(*p2, 'go', label='P2')
    plt.legend()
    plt.axis('equal')
    plt.title("Equal Arc Length Points on Quadratic Bézier Curve")
    plt.show()
