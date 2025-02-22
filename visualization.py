import numpy as np
import pandas as pd
pd.options.mode.chained_assignment = None
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.animation as animation

rng = np.random.default_rng(42)

TIME_STEPS = 160000
TIME_WRITE = 2000
NUM_FILES = int(TIME_STEPS / TIME_WRITE)
DEFECT_ORDER_THRESHOLD = 0.5  # lower is stricter
OBSTACLE_COLOR = 'xkcd:dark turquoise'

I = 420
J = 420
max_index = I * J - 1

h_chan_top_frac = 12 / 21
h_chan_bot_frac = 9 / 21
v_chan_left_frac = 9 / 21
v_chan_right_frac = 12 / 21

h_top_frac = 11 / 21
h_bot_frac = 10 / 21
h_left_frac = 3 / 21
h_right_frac = 18 / 21

v_top_frac = 18 / 21
v_bot_frac = 3 / 21
v_left_frac = 10 / 21
v_right_frac = 11 / 21

fig1, ax1 = plt.subplots()
fig1.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=None, hspace=None)
fig2, ax2 = plt.subplots()
fig2.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=None, hspace=None)

sample_fraction = 0.05
sample_indices = np.sort(rng.permutation(I * J)[:int(sample_fraction * I * J)])


def computeAngleDiff(startAngle, endAngle):
    min_index = np.argmin([np.abs(endAngle - startAngle), np.abs(endAngle - startAngle + np.pi), np.abs(endAngle - startAngle - np.pi)])

    if min_index == 0:
        return endAngle - startAngle
    elif min_index == 1:
        return endAngle - startAngle + np.pi
    else:
        return endAngle - startAngle - np.pi


def updateAx1(i):
    df_velocity = pd.read_csv('output/active_nematic_velocity_' + str(i * TIME_WRITE) + '.dat', sep=' ', header=None)
    x = df_velocity.iloc[:, 0]
    y = df_velocity.iloc[:, 1]
    vx = df_velocity.iloc[:, 3]
    vy = df_velocity.iloc[:, 4]

    ax1.clear()
    im = ax1.quiver(x[sample_indices], y[sample_indices], vx[sample_indices], vy[sample_indices], pivot='mid', scale=0.1, scale_units='width', width=0.0005, color='xkcd:royal blue')
    ax1.set_axis_off()
    xmin, xmax = np.min(x), np.max(x)
    ymin, ymax = np.min(y), np.max(y)

    # plus pattern
    h_point1 = np.array([xmin + h_left_frac * (xmax - xmin), ymin + h_bot_frac * (ymax - ymin)])
    h_point2 = np.array([xmin + h_left_frac * (xmax - xmin), ymin + h_top_frac * (ymax - ymin)])
    h_point3 = np.array([xmin + h_right_frac * (xmax - xmin), ymin + h_top_frac * (ymax - ymin)])
    h_point4 = np.array([xmin + h_right_frac * (xmax - xmin), ymin + h_bot_frac * (ymax - ymin)])

    v_point1 = np.array([xmin + v_left_frac * (xmax - xmin), ymin + v_bot_frac * (ymax - ymin)])
    v_point2 = np.array([xmin + v_left_frac * (xmax - xmin), ymin + v_top_frac * (ymax - ymin)])
    v_point3 = np.array([xmin + v_right_frac * (xmax - xmin), ymin + v_top_frac * (ymax - ymin)])
    v_point4 = np.array([xmin + v_right_frac * (xmax - xmin), ymin + v_bot_frac * (ymax - ymin)])

    corner_top_left = np.array([xmin + v_left_frac * (xmax - xmin), ymin + h_top_frac * (ymax - ymin)])
    corner_top_right = np.array([xmin + v_right_frac * (xmax - xmin), ymin + h_top_frac * (ymax - ymin)])
    corner_bot_right = np.array([xmin + v_right_frac * (xmax - xmin), ymin + h_bot_frac * (ymax - ymin)])
    corner_bot_left = np.array([xmin + v_left_frac * (xmax - xmin), ymin + h_bot_frac * (ymax - ymin)])

    plus = np.vstack((h_point2, corner_top_left, v_point2, v_point3, corner_top_right, h_point3, h_point4, corner_bot_right, v_point4, v_point1, corner_bot_left, h_point1))
    #channel = np.vstack((h_point2, h_point3, h_point4, h_point1))
    activity_pattern = patches.Polygon(plus, facecolor='xkcd:gold', alpha=0.3)
    ax1.add_patch(activity_pattern)

    rect_top_left = patches.Polygon(np.vstack((np.array([xmin, ymin + h_chan_top_frac * (ymax - ymin)]), np.array([xmin + v_chan_left_frac * (xmax - xmin), ymin + h_chan_top_frac * (ymax - ymin)]), np.array([xmin + v_chan_left_frac * (xmax - xmin), ymax]), np.array([xmin, ymax]))), facecolor=OBSTACLE_COLOR, alpha=1)
    ax1.add_patch(rect_top_left)

    rect_top_right = patches.Polygon(np.vstack((np.array([xmin + v_chan_right_frac * (xmax - xmin), ymin + h_chan_top_frac * (ymax - ymin)]), np.array([xmin + v_chan_right_frac * (xmax - xmin), ymax]), np.array([xmax, ymax]), np.array([xmax, ymin + h_chan_top_frac * (ymax - ymin)]))), facecolor=OBSTACLE_COLOR, alpha=1)
    ax1.add_patch(rect_top_right)

    rect_bot_right = patches.Polygon(np.vstack((np.array([xmax, ymin + h_chan_bot_frac * (ymax - ymin)]), np.array([xmax, ymin]), np.array([xmin + v_chan_right_frac * (xmax - xmin), ymin]), np.array([xmin + v_chan_right_frac * (xmax - xmin), ymin + h_chan_bot_frac * (ymax - ymin)]))), facecolor=OBSTACLE_COLOR, alpha=1)
    ax1.add_patch(rect_bot_right)

    rect_bot_left = patches.Polygon(np.vstack((np.array([xmin, ymin]), np.array([xmin, ymin + h_chan_bot_frac * (ymax - ymin)]), np.array([xmin + v_chan_left_frac * (xmax - xmin), ymin + h_chan_bot_frac * (ymax - ymin)]), np.array([xmin + v_chan_left_frac * (xmax - xmin), ymin]))), facecolor=OBSTACLE_COLOR, alpha=1)
    ax1.add_patch(rect_bot_left)

    ax1.axis('equal')

    return im,


def updateAx2(i):
    df_orientation = pd.read_csv('output/active_nematic_orientation_' + str(i * TIME_WRITE) + '.dat', sep=' ', header=None)
    x = df_orientation.iloc[:, 0]
    y = df_orientation.iloc[:, 1]
    order = df_orientation.iloc[:, 2]
    angle = df_orientation.iloc[:, 3]
    angle[angle < 0] += np.pi  # director goes from 0 to np.pi, not 0 to 2 * np.pi
    angle[angle > np.pi] -= np.pi
    cos = np.cos(angle)
    sin = np.sin(angle)

    x_defects = x.loc[order <= DEFECT_ORDER_THRESHOLD]
    y_defects = y.loc[order <= DEFECT_ORDER_THRESHOLD]

    filtered_x_defects = []
    filtered_y_defects = []
    defect_colors = []

    for x_defect, y_defect in zip(x_defects.to_numpy(), y_defects.to_numpy()):
        i0 = x_defect + y_defect * I
        i1 = x_defect - 1 + (y_defect + 1) * I
        i2 = x_defect + (y_defect + 1) * I
        i3 = x_defect + 1 + (y_defect + 1) * I
        i4 = x_defect - 1 + y_defect * I
        i5 = x_defect + 1 + y_defect * I
        i6 = x_defect - 1 + (y_defect - 1) * I
        i7 = x_defect + (y_defect - 1) * I
        i8 = x_defect + 1 + (y_defect - 1) * I
        defect_order = order[i0]

        if all(np.array([i1, i2, i3, i4, i5, i6, i7, i8]) < max_index) and all(np.array([i1, i2, i3, i4, i5, i6, i7, i8]) > 0):
            if defect_order < order[i1] and defect_order < order[i2] and defect_order < order[i3] and defect_order < order[i4] and defect_order < order[i5] and defect_order < order[i6] and defect_order < order[i7] and defect_order < order[i8]:
                winding_number = 0
                winding_number += computeAngleDiff(angle[i5], angle[i3])
                winding_number += computeAngleDiff(angle[i3], angle[i2])
                winding_number += computeAngleDiff(angle[i2], angle[i1])
                winding_number += computeAngleDiff(angle[i1], angle[i4])
                winding_number += computeAngleDiff(angle[i4], angle[i6])
                winding_number += computeAngleDiff(angle[i6], angle[i7])
                winding_number += computeAngleDiff(angle[i7], angle[i8])
                winding_number += computeAngleDiff(angle[i8], angle[i5])
                winding_number = np.round(winding_number / (2 * np.pi), decimals=1)
                
                if winding_number == 0.5:
                    filtered_x_defects.append(x_defect)
                    filtered_y_defects.append(y_defect)
                    defect_colors.append('xkcd:red')
                elif winding_number == -0.5:
                    filtered_x_defects.append(x_defect)
                    filtered_y_defects.append(y_defect)
                    defect_colors.append('xkcd:azure')

    ax2.clear()
    im = ax2.scatter(filtered_x_defects, filtered_y_defects, s=10, c=defect_colors)
    im = ax2.quiver(x[sample_indices], y[sample_indices], cos[sample_indices], sin[sample_indices], pivot='mid', width=0.0005, color='xkcd:royal blue', headlength=0, headaxislength=0)  # headless quivers for nematics
    ax2.set_axis_off()
    xmin, xmax = np.min(x), np.max(x)
    ymin, ymax = np.min(y), np.max(y)

    # plus pattern
    h_point1 = np.array([xmin + h_left_frac * (xmax - xmin), ymin + h_bot_frac * (ymax - ymin)])
    h_point2 = np.array([xmin + h_left_frac * (xmax - xmin), ymin + h_top_frac * (ymax - ymin)])
    h_point3 = np.array([xmin + h_right_frac * (xmax - xmin), ymin + h_top_frac * (ymax - ymin)])
    h_point4 = np.array([xmin + h_right_frac * (xmax - xmin), ymin + h_bot_frac * (ymax - ymin)])

    v_point1 = np.array([xmin + v_left_frac * (xmax - xmin), ymin + v_bot_frac * (ymax - ymin)])
    v_point2 = np.array([xmin + v_left_frac * (xmax - xmin), ymin + v_top_frac * (ymax - ymin)])
    v_point3 = np.array([xmin + v_right_frac * (xmax - xmin), ymin + v_top_frac * (ymax - ymin)])
    v_point4 = np.array([xmin + v_right_frac * (xmax - xmin), ymin + v_bot_frac * (ymax - ymin)])

    corner_top_left = np.array([xmin + v_left_frac * (xmax - xmin), ymin + h_top_frac * (ymax - ymin)])
    corner_top_right = np.array([xmin + v_right_frac * (xmax - xmin), ymin + h_top_frac * (ymax - ymin)])
    corner_bot_right = np.array([xmin + v_right_frac * (xmax - xmin), ymin + h_bot_frac * (ymax - ymin)])
    corner_bot_left = np.array([xmin + v_left_frac * (xmax - xmin), ymin + h_bot_frac * (ymax - ymin)])

    plus = np.vstack((h_point2, corner_top_left, v_point2, v_point3, corner_top_right, h_point3, h_point4, corner_bot_right, v_point4, v_point1, corner_bot_left, h_point1))
    #channel = np.vstack((h_point2, h_point3, h_point4, h_point1))
    activity_pattern = patches.Polygon(plus, facecolor='xkcd:gold', alpha=0.3)
    ax2.add_patch(activity_pattern)

    rect_top_left = patches.Polygon(np.vstack((np.array([xmin, ymin + h_chan_top_frac * (ymax - ymin)]), np.array([xmin + v_chan_left_frac * (xmax - xmin), ymin + h_chan_top_frac * (ymax - ymin)]), np.array([xmin + v_chan_left_frac * (xmax - xmin), ymax]), np.array([xmin, ymax]))), facecolor=OBSTACLE_COLOR, alpha=1)
    ax2.add_patch(rect_top_left)

    rect_top_right = patches.Polygon(np.vstack((np.array([xmin + v_chan_right_frac * (xmax - xmin), ymin + h_chan_top_frac * (ymax - ymin)]), np.array([xmin + v_chan_right_frac * (xmax - xmin), ymax]), np.array([xmax, ymax]), np.array([xmax, ymin + h_chan_top_frac * (ymax - ymin)]))), facecolor=OBSTACLE_COLOR, alpha=1)
    ax2.add_patch(rect_top_right)

    rect_bot_right = patches.Polygon(np.vstack((np.array([xmax, ymin + h_chan_bot_frac * (ymax - ymin)]), np.array([xmax, ymin]), np.array([xmin + v_chan_right_frac * (xmax - xmin), ymin]), np.array([xmin + v_chan_right_frac * (xmax - xmin), ymin + h_chan_bot_frac * (ymax - ymin)]))), facecolor=OBSTACLE_COLOR, alpha=1)
    ax2.add_patch(rect_bot_right)

    rect_bot_left = patches.Polygon(np.vstack((np.array([xmin, ymin]), np.array([xmin, ymin + h_chan_bot_frac * (ymax - ymin)]), np.array([xmin + v_chan_left_frac * (xmax - xmin), ymin + h_chan_bot_frac * (ymax - ymin)]), np.array([xmin + v_chan_left_frac * (xmax - xmin), ymin]))), facecolor=OBSTACLE_COLOR, alpha=1)
    ax2.add_patch(rect_bot_left)

    ax2.axis('equal')

    return im,


# Create the animation object
velocity_animation_fig = animation.FuncAnimation(fig1, updateAx1, frames=NUM_FILES, interval=40, blit=True, repeat_delay=2,)
orientation_animation_fig = animation.FuncAnimation(fig2, updateAx2, frames=NUM_FILES, interval=40, blit=True, repeat_delay=2,)
velocity_animation_fig.save("velocity.gif", dpi=400, savefig_kwargs=dict(facecolor='xkcd:white'))
orientation_animation_fig.save("orientation.gif", dpi=400, savefig_kwargs=dict(facecolor='xkcd:white'))