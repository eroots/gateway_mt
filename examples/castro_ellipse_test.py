"""Sketch of a rotated phase-tensor-style ellipse with axes and arrows.

Scratch figure script; does not use ``pyMT``.
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Polygon

# Rotation angle
angle_deg = 45
angle_rad = np.deg2rad(angle_deg)

# Center of the plot
center_x, center_y = 10, 5

# Rotation function
def rotate(x, y, angle):
    x_rot = x * np.cos(angle) - y * np.sin(angle)
    y_rot = x * np.sin(angle) + y * np.cos(angle)
    return x_rot, y_rot

# Ellipse parameters
a_large = 6
b_large = 4

# Arrow path (smaller ellipse)
scale = 0.7
a_arrow = a_large * scale
b_arrow = b_large * scale

# Full ellipse for plotting
theta_full = np.linspace(0, 2 * np.pi, 300)
x_large = a_large * np.cos(theta_full)
y_large = b_large * np.sin(theta_full)

# Partial arc for arrow: 20° to 45°
theta_arc = np.deg2rad(np.linspace(20, 45, 100))
x_arc = a_arrow * np.cos(theta_arc)
y_arc = b_arrow * np.sin(theta_arc)

# Rotate both shapes
x_large_rot, y_large_rot = rotate(x_large, y_large, angle_rad)
x_arc_rot, y_arc_rot = rotate(x_arc, y_arc, angle_rad)

# Translate to center (10, 5)
x_large_rot += center_x
y_large_rot += center_y
x_arc_rot += center_x
y_arc_rot += center_y

# Create filled arc patch
x_fill = np.concatenate(([center_x], x_arc_rot, [center_x]))
y_fill = np.concatenate(([center_y], y_arc_rot, [center_y]))
arc_patch = Polygon(np.column_stack((x_fill, y_fill)), closed=True, facecolor='red', edgecolor='none')

# Plotting
fig, ax = plt.subplots(figsize=(8, 6))

# Draw rotated ellipse
ax.plot(x_large_rot, y_large_rot, color='blue', label='Rotated Ellipse')

# Add filled arc
ax.add_patch(arc_patch)

# Draw arrow along the arc
arrow = FancyArrowPatch(
    (x_arc_rot[0], y_arc_rot[0]), (x_arc_rot[-1], y_arc_rot[-1]),
    connectionstyle="arc3,rad=0.3",
    arrowstyle="<|-", mutation_scale=20, color='red'
)
ax.add_patch(arrow)

# Formatting
ax.set_aspect('equal')
ax.set_xlim(center_x - a_large - 2, center_x + a_large + 2)
ax.set_ylim(center_y - b_large - 2, center_y + b_large + 2)
ax.set_title("Rotated Ellipse Centered at (10, 5) with Filled Partial Arc Arrow")
ax.axis('off')

plt.show()