import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.animation as animation

TIME_START = 100000
TIME_STEPS = 200000
TIME_WRITE = 400
NUM_FILES = int((TIME_STEPS - TIME_START) / TIME_WRITE) + 1
DEFECT_ORDER_THRESHOLD = 0.1  # lower is stricter

in_top_frac = 0.55
in_bot_frac = 0.45
in_left_frac = 0.05
in_right_frac = 0.95

out_top_frac = 0.45 #1
out_bot_frac = 0 #0.55
out_left_frac = 0.85
out_right_frac = 0.95

fig1, ax1 = plt.subplots()
fig1.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=None, hspace=None)
fig2, ax2 = plt.subplots()
fig2.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=None, hspace=None)


def updateAx1(i):
    df_velocity = pd.read_csv('output/active_nematic_velocity_' + str(TIME_START + i * TIME_WRITE) + '.dat', sep=' ', header=None)
    x = df_velocity.iloc[:, 0]
    y = df_velocity.iloc[:, 1]
    vx = df_velocity.iloc[:, 3]
    vy = df_velocity.iloc[:, 4]

    ax1.clear()
    im = ax1.quiver(x, y, vx, vy, pivot='mid', width=0.0005, color='xkcd:royal blue')
    ax1.set_axis_off()
    xmin, xmax = np.min(x), np.max(x)
    ymin, ymax = np.min(y), np.max(y)

    # inlet rectangle
    in_point1 = np.array([xmin + in_left_frac * (xmax - xmin), ymin + in_bot_frac * (ymax - ymin)])
    in_point2 = np.array([xmin + in_left_frac * (xmax - xmin), ymin + in_top_frac * (ymax - ymin)])
    in_point3 = np.array([xmin + in_right_frac * (xmax - xmin), ymin + in_top_frac * (ymax - ymin)])
    in_point4 = np.array([xmin + in_right_frac * (xmax - xmin), ymin + in_bot_frac * (ymax - ymin)])
    rectangle = np.vstack((in_point1, in_point2, in_point3, in_point4))
    activity_pattern = patches.Polygon(rectangle, facecolor='xkcd:gold', alpha=0.3)
    ax1.add_patch(activity_pattern)

    # outlet rectangle
    out_point1 = np.array([xmin + out_left_frac * (xmax - xmin), ymin + out_bot_frac * (ymax - ymin)])
    out_point2 = np.array([xmin + out_left_frac * (xmax - xmin), ymin + out_top_frac * (ymax - ymin)])
    out_point3 = np.array([xmin + out_right_frac * (xmax - xmin), ymin + out_top_frac * (ymax - ymin)])
    out_point4 = np.array([xmin + out_right_frac * (xmax - xmin), ymin + out_bot_frac * (ymax - ymin)])
    rectangle = np.vstack((out_point1, out_point2, out_point3, out_point4))
    activity_pattern = patches.Polygon(rectangle, facecolor='xkcd:gold', alpha=0.3)
    ax1.add_patch(activity_pattern)

    ax1.axis('equal')

    return im,

def updateAx2(i):
    df_orientation = pd.read_csv('output/active_nematic_orientation_' + str(TIME_START + i * TIME_WRITE) + '.dat', sep=' ', header=None)
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

    # inlet rectangle
    in_point1 = np.array([xmin + in_left_frac * (xmax - xmin), ymin + in_bot_frac * (ymax - ymin)])
    in_point2 = np.array([xmin + in_left_frac * (xmax - xmin), ymin + in_top_frac * (ymax - ymin)])
    in_point3 = np.array([xmin + in_right_frac * (xmax - xmin), ymin + in_top_frac * (ymax - ymin)])
    in_point4 = np.array([xmin + in_right_frac * (xmax - xmin), ymin + in_bot_frac * (ymax - ymin)])
    rectangle = np.vstack((in_point1, in_point2, in_point3, in_point4))
    activity_pattern = patches.Polygon(rectangle, facecolor='xkcd:gold', alpha=0.3)
    ax2.add_patch(activity_pattern)

    # outlet rectangle
    out_point1 = np.array([xmin + out_left_frac * (xmax - xmin), ymin + out_bot_frac * (ymax - ymin)])
    out_point2 = np.array([xmin + out_left_frac * (xmax - xmin), ymin + out_top_frac * (ymax - ymin)])
    out_point3 = np.array([xmin + out_right_frac * (xmax - xmin), ymin + out_top_frac * (ymax - ymin)])
    out_point4 = np.array([xmin + out_right_frac * (xmax - xmin), ymin + out_bot_frac * (ymax - ymin)])
    rectangle = np.vstack((out_point1, out_point2, out_point3, out_point4))
    activity_pattern = patches.Polygon(rectangle, facecolor='xkcd:gold', alpha=0.3)
    ax2.add_patch(activity_pattern)

    ax2.axis('equal')

    return im,

# Create the animation object
velocity_animation_fig = animation.FuncAnimation(fig1, updateAx1, frames=NUM_FILES, interval=20, blit=True, repeat_delay=2,)
orientation_animation_fig = animation.FuncAnimation(fig2, updateAx2, frames=NUM_FILES, interval=20, blit=True, repeat_delay=2,)
velocity_animation_fig.save("velocity.gif", dpi=400, savefig_kwargs=dict(facecolor='xkcd:white'))
orientation_animation_fig.save("orientation.gif", dpi=400, savefig_kwargs=dict(facecolor='xkcd:white'))