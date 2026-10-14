"""Cell-average initialization and an independent smooth exact benchmark."""

import numpy as np


def cell_average(grid, function, order=8):
    """Gauss-Legendre cell averages, not point samples at cell centers."""
    nodes, weights = np.polynomial.legendre.leggauss(order)
    points = grid.x[:, None] + 0.5 * grid.dx * nodes[None, :]
    values = np.asarray(function(points), dtype=float)
    if values.shape[:2] != points.shape:
        raise ValueError("Function must preserve the quadrature input shape")
    return 0.5 * np.einsum("j,ij...->i...", weights, values)


def right_going_state(h, depth, g=9.81):
    """Constant-bed nonlinear simple wave: u-2 sqrt(gh)=-2 sqrt(gD)."""
    h = np.asarray(h, dtype=float)
    u = 2 * (np.sqrt(g * h) - np.sqrt(g * depth))
    return np.stack((h, h * u), axis=-1)


def sine_simple_wave(x, time, depth, amplitude, length, g=9.81):
    """Exact smooth periodic nonlinear simple wave, before characteristic crossing.

    Solve x = xi + (3 sqrt(g h0(xi))-2 sqrt(gD)) t using Newton.
    This uses characteristics, not the finite-volume solver being tested.
    A conservative sufficient no-crossing condition is checked explicitly.
    """
    if not (0 <= amplitude < depth and length > 0 and depth > 0 and g > 0 and time >= 0):
        raise ValueError("Invalid exact-wave parameters")
    k, c0 = 2 * np.pi / length, np.sqrt(g * depth)
    worst_slope = 1.5 * np.sqrt(g / (depth - amplitude)) * amplitude * k
    if time * worst_slope >= 0.95:
        raise ValueError("Exact characteristic benchmark too close to wave breaking")
    x = np.asarray(x, dtype=float)
    foot = x - c0 * time
    for _ in range(30):
        h = depth + amplitude * np.sin(k * foot)
        c = np.sqrt(g * h)
        speed = 3 * c - 2 * c0
        derivative = 1 + time * 1.5 * g * amplitude * k * np.cos(k * foot) / c
        update = (foot + speed * time - x) / derivative
        foot -= update
        if np.max(np.abs(update)) < 1e-12 * max(1.0, length):
            break
    else:
        raise RuntimeError("Exact characteristic solve did not converge")
    h = depth + amplitude * np.sin(k * foot)
    return right_going_state(h, depth, g)


def wet_dam_break(x, time, left_depth=2.0, right_depth=1.0, g=9.81):
    """Independent exact wet-bed dam break: left rarefaction, right shock.

    Initial velocity is zero and the interface is at x=0. This benchmark is
    not a tsunami initial condition. A scalar bisection gives the star state.
    """
    if not (left_depth > right_depth > 0 and g > 0 and time >= 0):
        raise ValueError("Require left_depth > right_depth > 0 and time >= 0")
    x = np.asarray(x, dtype=float)
    if time == 0:
        h = np.where(x < 0, left_depth, right_depth)
        return np.stack((h, np.zeros_like(h)), axis=-1)
    cl = np.sqrt(g * left_depth)
    lower, upper = right_depth, left_depth
    for _ in range(80):
        hm = 0.5 * (lower + upper)
        mismatch = 2*(np.sqrt(g*hm)-cl) + (hm-right_depth)*np.sqrt(0.5*g*(1/hm+1/right_depth))
        if mismatch > 0:
            upper = hm
        else:
            lower = hm
    hm = 0.5*(lower+upper)
    um = 2*(cl-np.sqrt(g*hm))
    shock_speed = hm*um/(hm-right_depth)
    s = x/time
    h = np.full_like(s, right_depth)
    u = np.zeros_like(s)
    left = s <= -cl
    fan = (s > -cl) & (s < um-np.sqrt(g*hm))
    middle = (s >= um-np.sqrt(g*hm)) & (s < shock_speed)
    h[left] = left_depth
    h[fan] = ((2*cl-s[fan])/3)**2/g
    u[fan] = 2*(cl+s[fan])/3
    h[middle], u[middle] = hm, um
    return np.stack((h, h*u), axis=-1)
