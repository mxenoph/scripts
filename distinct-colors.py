#!/usr/bin/env python

import argparse
import importlib

parser = argparse.ArgumentParser(description='Find suitable palette\n')
parser.add_argument('-n', '--number', required = True, help = 'Number of colors to plot')
parser.add_argument('-o', '--out_dir', required = True, help = 'Directory to save the plot')
args = parser.parse_args()

def find_palette(n, palette = 'qualitative', sets = 'Set3_'):
    if palette == 'diverging':
        sets = 'RdBu_'
    elif palette == 'sequential':
        if int(n) > 9 :
            palette = 'diverging'
            sets = 'RdBu_'
        else:
            sets = 'YlGnBu_'

    get = 'palettable.colorbrewer.' + palette
    colors = sets + str(n)
    palettes = importlib.import_module(get)
    p = getattr(palettes, colors)
    return p

def sample_plot(color_cycle, n, pp):
    import numpy as np
    import matplotlib.pyplot as plt
    # evenly sampled time at 200ms intervals
    t = np.arange(0., 5., 0.2)

    fig = plt.figure()
    ax = fig.add_subplot(111)

    fig.gca().set_color_cycle(a.mpl_colors)
#    plt.set_color_cycle(color_cycle)
    # red dashes, blue squares and green triangles
    ax.plot(t, t, t, t** 2, t, t ** 3, t, t**4, t, t**5)
    pp.savefig(fig)
    return plt

from matplotlib.backends.backend_pdf import PdfPages
pp = PdfPages(args.out_dir + '/find_palette.pdf')
a = find_palette(n = args.number, palette = 'sequential')
sample_plot(a, args.number, pp)
pp.close()
