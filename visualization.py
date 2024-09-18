import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.animation as animation

TIME_STEPS = 30000
TIME_WRITE = 50
NUM_FILES = int(TIME_STEPS / TIME_WRITE)
DEFECT_ORDER_THRESHOLD = 0.7
ACTIVITY_X_WIDTH = 0.9
ACTIVITY_Y_WIDTH = 0.2
ROTATION_ANGLE = 0 #0.25 * np.pi

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
    im = ax1.quiver(x, y, vx, vy, width=0.0005, color='xkcd:red orange')
    ax1.set_axis_off()
    xmin, xmax = np.min(x), np.max(x)
    ymin, ymax = np.min(y), np.max(y)
    x_anchor = (xmin + xmax) / 2 - ACTIVITY_X_WIDTH * (xmax - xmin) / 2
    y_anchor = (ymin + ymax) / 2 - ACTIVITY_Y_WIDTH * (ymax - ymin) / 2
    rotation_point = (xmin + 0.2 * (xmax - xmin), (ymin + ymax) / 2)
    activity_pattern = patches.Rectangle((x_anchor, y_anchor), ACTIVITY_X_WIDTH * (xmax - xmin), ACTIVITY_Y_WIDTH * (ymax - ymin), rotation_point=rotation_point, angle=np.rad2deg(ROTATION_ANGLE), facecolor='xkcd:gold', alpha=0.3)
    ax1.add_patch(activity_pattern)
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

    ax2.clear()
    im = ax2.scatter(x_defects, y_defects, c='xkcd:azure')
    im = ax2.quiver(x, y, cos, sin, width=0.0005, color='xkcd:red orange', headlength=0, headaxislength=0)  # headless quivers for nematics
    ax2.set_axis_off()
    xmin, xmax = np.min(x), np.max(x)
    ymin, ymax = np.min(y), np.max(y)
    x_anchor = (xmin + xmax) / 2 - ACTIVITY_X_WIDTH * (xmax - xmin) / 2
    y_anchor = (ymin + ymax) / 2 - ACTIVITY_Y_WIDTH * (ymax - ymin) / 2
    rotation_point = (xmin + 0.2 * (xmax - xmin), (ymin + ymax) / 2)
    activity_pattern = patches.Rectangle((x_anchor, y_anchor), ACTIVITY_X_WIDTH * (xmax - xmin), ACTIVITY_Y_WIDTH * (ymax - ymin), rotation_point=rotation_point, angle=np.rad2deg(ROTATION_ANGLE), facecolor='xkcd:gold', alpha=0.3)
    ax2.add_patch(activity_pattern)
    return im,

# Create the animation object
velocity_animation_fig = animation.FuncAnimation(fig1, updateAx1, frames=NUM_FILES, interval=100, blit=True, repeat_delay=2,)
orientation_animation_fig = animation.FuncAnimation(fig2, updateAx2, frames=NUM_FILES, interval=100, blit=True, repeat_delay=2,)
velocity_animation_fig.save("velocity.gif", dpi=300, savefig_kwargs=dict(facecolor='xkcd:almost black'))
orientation_animation_fig.save("orientation.gif", dpi=300, savefig_kwargs=dict(facecolor='xkcd:almost black'))