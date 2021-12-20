import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Arc
from matplotlib import cm
from matplotlib.colors import rgb2hex, Normalize
import matplotlib
import pandas as pd
import numpy as np

matplotlib.rcParams.update({'font.size': 14})


def draw_model(n_passes=8, n_elements=7, show=False, image_file=None):

    fig, ax = plt.subplots(figsize=(9, 7))

    wd = 1000
    ht = 1000
    lw = 2

    edge_pad = 100
    mid_pad = 40

    ax.set_xlim(0, wd)
    ax.set_ylim(0, ht)

    ele_wd = (wd - edge_pad * 2 - mid_pad * (n_elements - 1)) / n_elements
    ele_ht = (ht - edge_pad * 2 - mid_pad * (n_passes - 1)) / n_passes

    # sCO2 arrow in...
    ax.arrow(0, ht - (edge_pad + ele_ht / 2), dx=edge_pad, dy=0,
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

    ax.arrow(arrow_x, edge_pad + ele_ht / 2, dx=arrow_dx, dy=0,
             length_includes_head=True,
             width=lw,
             head_width=lw*10,
             color='k')

    # Air side arrows...
    arrow_x = edge_pad + ele_wd / 2
    for i in range(n_elements):
        ax.arrow(arrow_x, ht-edge_pad, dx=0, dy=edge_pad/2,
                 length_includes_head=True,
                 width=lw,
                 head_width=lw * 10,
                 color='k')
        ax.arrow(arrow_x, edge_pad/2, dx=0, dy=edge_pad/2,
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

        ax.arrow(edge_pad, y + ele_ht / 2, dx=wd-edge_pad*2, dy=0, width=lw)

        for i in range(n_elements):

            tube_rect = Rectangle((x, y), ele_wd, ele_ht, edgecolor='k', facecolor='w')

            ax.add_patch(tube_rect)

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
            ax.add_patch(arc)

        x_inc = -x_inc
        x += x_inc
        y += mid_pad + ele_ht

    ax.axis('off')

    if image_file is not None:
        fig.savefig(image_file)
    if show:
        plt.show()


draw_model(8, 7, show=True)
