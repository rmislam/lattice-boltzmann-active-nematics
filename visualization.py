import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.animation as animation

TIME_STEPS = 20000
TIME_WRITE = 100
NUM_FILES = int(TIME_STEPS / TIME_WRITE)
DEFECT_ORDER_THRESHOLD = 0.1  # lower is stricter

tri1_base_x_frac = 0.1
tri2_base_x_frac = 0.27
tri3_base_x_frac = 0.44
tri4_base_x_frac = 0.61
tri5_base_x_frac = 0.78
tri_base_y_frac = 0.5
height_frac = 0.15
half_width_frac = 0.1

fig1, ax1 = plt.subplots()
fig1.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=None, hspace=None)
fig2, ax2 = plt.subplots()
fig2.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=None, hspace=None)


def updateAx1(i):
    df_velocity = pd.read_csv('output/active_nematic_velocity_' + str(i * TIME_WRITE) + '.dat', sep=' ', header=None)
    x = df_velocity.iloc[:, 0]
    y = df_velocity.iloc[:, 1]
    vx = df_velocity.iloc[:, 3]
    vy = df_velocity.iloc[:, 4]

    ax1.clear()
    im = ax1.quiver(x, y, vx, vy, pivot='mid', width=0.0005, color='xkcd:royal blue')
    ax1.set_axis_off()
    xmin, xmax = np.min(x), np.max(x)
    ymin, ymax = np.min(y), np.max(y)

    # triangle 1
    bottom_point = np.array([xmin + tri1_base_x_frac * (xmax - xmin), ymin + (tri_base_y_frac + half_width_frac) * (ymax - ymin)])  # bottom point
    top_point = np.array([xmin + tri1_base_x_frac * (xmax - xmin), ymin + (tri_base_y_frac - half_width_frac) * (ymax - ymin)])  # top point
    right_point = np.array([xmin + (tri1_base_x_frac + height_frac) * (xmax - xmin), ymin + tri_base_y_frac * (ymax - ymin)])  # right point
    triangle = np.vstack((bottom_point, top_point, right_point))
    activity_pattern = patches.Polygon(triangle, facecolor='xkcd:gold', alpha=0.3)
    ax1.add_patch(activity_pattern)

    # triangle 2
    bottom_point = np.array([xmin + tri2_base_x_frac * (xmax - xmin), ymin + (tri_base_y_frac + half_width_frac) * (ymax - ymin)])  # bottom point
    top_point = np.array([xmin + tri2_base_x_frac * (xmax - xmin), ymin + (tri_base_y_frac - half_width_frac) * (ymax - ymin)])  # top point
    right_point = np.array([xmin + (tri2_base_x_frac + height_frac) * (xmax - xmin), ymin + tri_base_y_frac * (ymax - ymin)])  # right point
    triangle = np.vstack((bottom_point, top_point, right_point))
    activity_pattern = patches.Polygon(triangle, facecolor='xkcd:gold', alpha=0.3)
    ax1.add_patch(activity_pattern)

    # triangle 3
    bottom_point = np.array([xmin + tri3_base_x_frac * (xmax - xmin), ymin + (tri_base_y_frac + half_width_frac) * (ymax - ymin)])  # bottom point
    top_point = np.array([xmin + tri3_base_x_frac * (xmax - xmin), ymin + (tri_base_y_frac - half_width_frac) * (ymax - ymin)])  # top point
    right_point = np.array([xmin + (tri3_base_x_frac + height_frac) * (xmax - xmin), ymin + tri_base_y_frac * (ymax - ymin)])  # right point
    triangle = np.vstack((bottom_point, top_point, right_point))
    activity_pattern = patches.Polygon(triangle, facecolor='xkcd:gold', alpha=0.3)
    ax1.add_patch(activity_pattern)

    # triangle 4
    bottom_point = np.array([xmin + tri4_base_x_frac * (xmax - xmin), ymin + (tri_base_y_frac + half_width_frac) * (ymax - ymin)])  # bottom point
    top_point = np.array([xmin + tri4_base_x_frac * (xmax - xmin), ymin + (tri_base_y_frac - half_width_frac) * (ymax - ymin)])  # top point
    right_point = np.array([xmin + (tri4_base_x_frac + height_frac) * (xmax - xmin), ymin + tri_base_y_frac * (ymax - ymin)])  # right point
    triangle = np.vstack((bottom_point, top_point, right_point))
    activity_pattern = patches.Polygon(triangle, facecolor='xkcd:gold', alpha=0.3)
    ax1.add_patch(activity_pattern)

    # triangle 5
    bottom_point = np.array([xmin + tri5_base_x_frac * (xmax - xmin), ymin + (tri_base_y_frac + half_width_frac) * (ymax - ymin)])  # bottom point
    top_point = np.array([xmin + tri5_base_x_frac * (xmax - xmin), ymin + (tri_base_y_frac - half_width_frac) * (ymax - ymin)])  # top point
    right_point = np.array([xmin + (tri5_base_x_frac + height_frac) * (xmax - xmin), ymin + tri_base_y_frac * (ymax - ymin)])  # right point
    triangle = np.vstack((bottom_point, top_point, right_point))
    activity_pattern = patches.Polygon(triangle, facecolor='xkcd:gold', alpha=0.3)
    ax1.add_patch(activity_pattern)

    ax1.axis('equal')

    return im,

def updateAx2(i):
    df_orientation = pd.read_csv('output/active_nematic_orientation_' + str(i * TIME_WRITE) + '.dat', sep=' ', header=None)
    x = df_orientation.iloc[:, 0]
    y = df_orientation.iloc[:, 1]
    order = df_orientation.iloc[:, 2]
    angle = df_orientation.iloc[:, 3]
    cos = np.cos(angle)
    sin = np.sin(angle)

    x_defects = x.loc[order <= DEFECT_ORDER_THRESHOLD]
    y_defects = y.loc[order <= DEFECT_ORDER_THRESHOLD]
    #from IPython import embed; embed()

    ax2.clear()
    im = ax2.scatter(x_defects, y_defects, c='xkcd:azure')
    im = ax2.quiver(x, y, cos, sin, pivot='mid', width=0.0005, color='xkcd:royal blue', headlength=0, headaxislength=0)  # headless quivers for nematics
    ax2.set_axis_off()
    xmin, xmax = np.min(x), np.max(x)
    ymin, ymax = np.min(y), np.max(y)

    # triangle 1
    bottom_point = np.array([xmin + tri1_base_x_frac * (xmax - xmin), ymin + (tri_base_y_frac + half_width_frac) * (ymax - ymin)])  # bottom point
    top_point = np.array([xmin + tri1_base_x_frac * (xmax - xmin), ymin + (tri_base_y_frac - half_width_frac) * (ymax - ymin)])  # top point
    right_point = np.array([xmin + (tri1_base_x_frac + height_frac) * (xmax - xmin), ymin + tri_base_y_frac * (ymax - ymin)])  # right point
    triangle = np.vstack((bottom_point, top_point, right_point))
    activity_pattern = patches.Polygon(triangle, facecolor='xkcd:gold', alpha=0.3)
    ax2.add_patch(activity_pattern)

    # triangle 2
    bottom_point = np.array([xmin + tri2_base_x_frac * (xmax - xmin), ymin + (tri_base_y_frac + half_width_frac) * (ymax - ymin)])  # bottom point
    top_point = np.array([xmin + tri2_base_x_frac * (xmax - xmin), ymin + (tri_base_y_frac - half_width_frac) * (ymax - ymin)])  # top point
    right_point = np.array([xmin + (tri2_base_x_frac + height_frac) * (xmax - xmin), ymin + tri_base_y_frac * (ymax - ymin)])  # right point
    triangle = np.vstack((bottom_point, top_point, right_point))
    activity_pattern = patches.Polygon(triangle, facecolor='xkcd:gold', alpha=0.3)
    ax2.add_patch(activity_pattern)

    # triangle 3
    bottom_point = np.array([xmin + tri3_base_x_frac * (xmax - xmin), ymin + (tri_base_y_frac + half_width_frac) * (ymax - ymin)])  # bottom point
    top_point = np.array([xmin + tri3_base_x_frac * (xmax - xmin), ymin + (tri_base_y_frac - half_width_frac) * (ymax - ymin)])  # top point
    right_point = np.array([xmin + (tri3_base_x_frac + height_frac) * (xmax - xmin), ymin + tri_base_y_frac * (ymax - ymin)])  # right point
    triangle = np.vstack((bottom_point, top_point, right_point))
    activity_pattern = patches.Polygon(triangle, facecolor='xkcd:gold', alpha=0.3)
    ax2.add_patch(activity_pattern)

    # triangle 4
    bottom_point = np.array([xmin + tri4_base_x_frac * (xmax - xmin), ymin + (tri_base_y_frac + half_width_frac) * (ymax - ymin)])  # bottom point
    top_point = np.array([xmin + tri4_base_x_frac * (xmax - xmin), ymin + (tri_base_y_frac - half_width_frac) * (ymax - ymin)])  # top point
    right_point = np.array([xmin + (tri4_base_x_frac + height_frac) * (xmax - xmin), ymin + tri_base_y_frac * (ymax - ymin)])  # right point
    triangle = np.vstack((bottom_point, top_point, right_point))
    activity_pattern = patches.Polygon(triangle, facecolor='xkcd:gold', alpha=0.3)
    ax2.add_patch(activity_pattern)

    # triangle 5
    bottom_point = np.array([xmin + tri5_base_x_frac * (xmax - xmin), ymin + (tri_base_y_frac + half_width_frac) * (ymax - ymin)])  # bottom point
    top_point = np.array([xmin + tri5_base_x_frac * (xmax - xmin), ymin + (tri_base_y_frac - half_width_frac) * (ymax - ymin)])  # top point
    right_point = np.array([xmin + (tri5_base_x_frac + height_frac) * (xmax - xmin), ymin + tri_base_y_frac * (ymax - ymin)])  # right point
    triangle = np.vstack((bottom_point, top_point, right_point))
    activity_pattern = patches.Polygon(triangle, facecolor='xkcd:gold', alpha=0.3)
    ax2.add_patch(activity_pattern)

    ax2.axis('equal')

    return im,

# Create the animation object
velocity_animation_fig = animation.FuncAnimation(fig1, updateAx1, frames=NUM_FILES, interval=20, blit=True, repeat_delay=2,)
orientation_animation_fig = animation.FuncAnimation(fig2, updateAx2, frames=NUM_FILES, interval=20, blit=True, repeat_delay=2,)
velocity_animation_fig.save("velocity.gif", dpi=400, savefig_kwargs=dict(facecolor='xkcd:white'))
orientation_animation_fig.save("orientation.gif", dpi=400, savefig_kwargs=dict(facecolor='xkcd:white'))