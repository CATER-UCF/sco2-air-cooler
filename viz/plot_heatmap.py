"""
Plots the two-dimensional temperature heatmap used in the paper. Note: the image
should be manually cropped to remove some extra whitespace.
"""
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Arc
from matplotlib import cm, rc
from matplotlib.colors import rgb2hex, Normalize
import pandas as pd
import numpy as np

rc('font', **{'family': 'sans-serif', 'sans-serif': ['DejaVu Sans'], 'size': 10})
rc('mathtext', **{'default': 'regular'})


def plot_2d_steady_state(data_file, n_passes=8, n_elements=7, show=False, image_file=None):

    # Read in the results data
    df = pd.read_csv(data_file)
    n_eles = n_passes * n_elements
    tube_temps, shell_temps = np.empty(n_eles), np.empty(n_eles)
    for i in range(n_eles):
        tube_series = df[f'temperature_tube_in_{i}']
        shell_series = df[f'temperature_shell_in_{i}']
        tube_temps[i] = tube_series[0]
        shell_temps[i] = shell_series[0]

    # Since CO2 flows from top to bottom and air vice versa...
    tube_temps = np.flip(tube_temps)
    shell_temps = np.flip(shell_temps)

    min_temp = min(np.min(shell_temps), min(tube_temps))
    max_temp = max(np.max(shell_temps), max(tube_temps))

    color_res = 1000
    cmap = cm.get_cmap('turbo', color_res)

    fig, ax = plt.subplots(1, 2, figsize=(10, 4))

    wd = 1000
    ht = 1000
    lw = 2

    edge_pad = 100
    mid_pad = 40

    ax[0].set_xlim(0, wd)
    ax[0].set_ylim(0, ht)
    ax[1].set_xlim(0, wd)
    ax[1].set_ylim(0, ht)

    ele_wd = (wd - edge_pad * 2 - mid_pad * (n_elements - 1)) / n_elements
    ele_ht = (ht - edge_pad * 2 - mid_pad * (n_passes - 1)) / n_passes

    # sCO2 arrow in...
    ax[0].arrow(0, ht - (edge_pad + ele_ht / 2), dx=edge_pad, dy=0,
                length_includes_head=True,
                width=lw,
                head_width=lw*10,
                color='k')

    # sCO2 arrow out...
    if n_passes % 2:
        arrow_x = wd - edge_pad
        arrow_dx = edge_pad
    else:
        arrow_x = edge_pad
        arrow_dx = -edge_pad

    ax[0].arrow(arrow_x, edge_pad + ele_ht / 2, dx=arrow_dx, dy=0,
                length_includes_head=True,
                width=lw,
                head_width=lw*10,
                color='k')

    # Air side arrows...
    arrow_x = edge_pad + ele_wd / 2
    for i in range(n_elements):
        ax[1].arrow(arrow_x, edge_pad, dx=0, dy=ht-3*edge_pad/2,
                    length_includes_head=True,
                    width=lw,
                    head_width=lw * 10,
                    color='k')
        ax[1].arrow(arrow_x, edge_pad/2, dx=0, dy=edge_pad/2,
                    length_includes_head=True,
                    width=lw,
                    head_width=lw * 10,
                    color='k')

        arrow_x += mid_pad + ele_wd

    # Coordinates of element 0
    x, y = edge_pad, edge_pad

    x_inc = mid_pad + ele_wd
    ele_idx = 0

    for j in range(n_passes):

        ax[0].arrow(edge_pad, y + ele_ht / 2, dx=wd-edge_pad*2, dy=0, width=lw)

        for i in range(n_elements):

            tube_color = cmap(int((color_res * (tube_temps[ele_idx] - min_temp)) / (max_temp - min_temp)))
            tube_rect = Rectangle((x, y), ele_wd, ele_ht, color=rgb2hex(tube_color))

            shell_color = cmap(int((color_res * (shell_temps[ele_idx] - min_temp)) / (max_temp - min_temp)))
            shell_rect = Rectangle((x, y), ele_wd, ele_ht, color=rgb2hex(shell_color))

            ax[0].add_patch(tube_rect)
            ax[1].add_patch(shell_rect)

            x += x_inc
            ele_idx += 1

        # At each pass, add a turn...
        arc_y = y + (mid_pad / 2 + ele_ht)
        diameter = mid_pad + ele_ht
        if x_inc > 0:
            arc_x = wd - edge_pad
            th1 = 270
            th2 = 90
        else:
            arc_x = edge_pad
            th1 = 90
            th2 = 270

        if j < n_passes - 1:
            arc = Arc((arc_x, arc_y), width=1.2*diameter, height=diameter, theta1=th1, theta2=th2, linewidth=lw)
            ax[0].add_patch(arc)

        x_inc = -x_inc
        x += x_inc
        y += mid_pad + ele_ht

    ax[0].axis('off')
    ax[1].axis('off')

    ax[0].set_title('sCO2 Flow')
    ax[1].set_title('Air Flow')

    anorm = Normalize(vmin=min_temp-273.15, vmax=max_temp-273.15)
    cb = fig.colorbar(cm.ScalarMappable(norm=anorm, cmap=cmap), ax=ax, shrink=0.8)
    cb.set_label('Temperature (C)', rotation=270, labelpad=20)

    if image_file is not None:
        fig.savefig(image_file, dpi=500)
    if show:
        plt.show()


plot_2d_steady_state('./data/combined.csv',
                     image_file='./images/temperature_profile_2d.png', show=True)
