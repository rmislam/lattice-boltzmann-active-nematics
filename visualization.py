import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation

TIME_STEPS = 10000
TIME_WRITE = 50
NUM_FILES = int(TIME_STEPS / TIME_WRITE)
DEFECT_ORDER_THRESHOLD = 0.7

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
    im = ax1.quiver(x, y, vx, vy, width=0.001, color='blue')
    ax1.set_axis_off()
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
    im = ax2.scatter(x_defects, y_defects, c='green')
    im = ax2.quiver(x, y, cos, sin, width=0.001, color='blue')
    ax2.set_axis_off()
    return im,

# Create the animation object
velocity_animation_fig = animation.FuncAnimation(fig1, updateAx1, frames=NUM_FILES, interval=100, blit=True, repeat_delay=5,)
orientation_animation_fig = animation.FuncAnimation(fig2, updateAx2, frames=NUM_FILES, interval=100, blit=True, repeat_delay=5,)
velocity_animation_fig.save("velocity.gif", dpi=600)
orientation_animation_fig.save("orientation.gif", dpi=600)