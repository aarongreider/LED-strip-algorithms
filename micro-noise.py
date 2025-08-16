from random import uniform
import math
import plasma
import time

# ----------------------------
# CONFIG
# ----------------------------
NUM_LEDS = 144
FADE_UP_SPEED = 2
FADE_DOWN_SPEED = 2

# ----------------------------
# BEZIER + NOISE FUNCTIONS
# ----------------------------
def quad_bezier_point(p0, p1, p2, t):
    u = 1 - t
    return (
        u*u*p0[0] + 2*u*t*p1[0] + t*t*p2[0],
        u*u*p0[1] + 2*u*t*p1[1] + t*t*p2[1]
    )

def approx_bezier_length(p0, p1, p2, steps=20):
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
    lengths = []
    total_length = 0
    for p0, p1, p2 in curves:
        seg_len = approx_bezier_length(p0, p1, p2, steps=steps_per_curve)
        lengths.append(seg_len)
        total_length += seg_len

    spacing = total_length / (num_points - 1)
    points = []
    seg_index = 0
    dist_in_seg = 0
    p0, p1, p2 = curves[0]
    prev_point = quad_bezier_point(p0, p1, p2, 0)
    points.append(prev_point)
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
            if dist_in_seg >= spacing:
                points.append(curr_point)
                dist_in_seg = 0
                prev_point = curr_point
                if len(points) >= num_points:
                    break
            else:
                prev_point = curr_point
        seg_index += 1
        t_index = 0
        dist_in_seg = 0
        if seg_index < len(curves):
            prev_point = quad_bezier_point(*curves[seg_index], 0)
    while len(points) < num_points:
        points.append(points[-1])
    return points

def random2(st):
    dot = st[0] * 12.9898 + st[1] * 78.233
    return (math.sin(dot) * 43758.5453123) % 1.0

def noise2(st):
    i = [math.floor(st[0]), math.floor(st[1])]
    f = [st[0] - i[0], st[1] - i[1]]
    a = random2(i)
    b = random2([i[0] + 1.0, i[1]])
    c = random2([i[0], i[1] + 1.0])
    d = random2([i[0] + 1.0, i[1] + 1.0])
    u = [f[0]*f[0]*(3.0 - 2.0*f[0]), f[1]*f[1]*(3.0 - 2.0*f[1])]
    n = (a * (1 - u[0]) + b * u[0]) * (1 - u[1]) + (c * (1 - u[0]) + d * u[0]) * u[1]
    return n

# ----------------------------
# LED HANDLING
# ----------------------------
current_leds = [[0] * 3 for _ in range(NUM_LEDS)]
target_leds = [[0] * 3 for _ in range(NUM_LEDS)]

led_strip = plasma.WS2812(NUM_LEDS, color_order=plasma.COLOR_ORDER_RGB)
led_strip.start()

def display_current():
    for i in range(NUM_LEDS):
        r, g, b = current_leds[i]
        led_strip.set_rgb(i, int(r), int(g), int(b))

def move_to_target():
    for i in range(NUM_LEDS):
        for c in range(3):
            if current_leds[i][c] < target_leds[i][c]:
                current_leds[i][c] = min(current_leds[i][c] + FADE_UP_SPEED, target_leds[i][c])
            elif current_leds[i][c] > target_leds[i][c]:
                current_leds[i][c] = max(current_leds[i][c] - FADE_DOWN_SPEED, target_leds[i][c])

# ----------------------------
# MAIN
# ----------------------------
curves = [
    ((0, 0), (10, 0), (10, 10)),
    ((10, 10), (6, 14), (2, 10)),
    ((2, 10), (2, 0), (12, 0)),
]
points = stitch_quadratic_beziers(curves, NUM_LEDS, steps_per_curve=200)

frame = 0
last_time = time.ticks_ms()  # record the first timestamp

while True:
    print(frame)
    t = frame * 0.05
    for i, p in enumerate(points):
        n = noise2([p[0] + t, p[1] + t * 0.3])
        r = int(n * 255)
        g = int((1 - n) * 255)
        b = 128
        target_leds[i] = [r, g, b]
    move_to_target()
    display_current()
    frame += 1
    #time.sleep(0.05)  # ~20 FPS
