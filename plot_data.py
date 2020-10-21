#!/usr/bin/env python3
"""
Plots the results
"""

import numpy as np
import matplotlib.pyplot as plt
from itertools import islice


def read_data():
    with open('route.dat', 'r') as f:
        n = f.readline()
        n = int(n.strip('\n'))
        rawdata = [x.strip('\n').split(' ') for x in islice(f, int(n))]

        coordinates = np.array(
            list(
                zip(np.array([x for x, _ in rawdata]), np.array([y for _, y in rawdata]))
            ), dtype=[('x', np.float), ('y', np.float)])

    rdata = np.loadtxt("data/racetrap.data")

    return rdata, coordinates


def plot_data():
    rdata, coordinates = read_data()
    dim = rdata.shape
    i = 0
    for path in rdata:
        i = i + 1
        path = path.tolist()
        path.append(0)

        plt.plot([coordinates[int(x)]['x'] for x in path], [coordinates[int(y)]['y'] for y in path], 'r', zorder=1,
                 lw=3)
        plt.scatter([coordinates[int(x)]['x'] for x in path], [coordinates[int(y)]['y'] for y in path], zorder=2)
        plt.plot()
        plt.pause(0.005)

        if i == dim[0]:
            plt.pause(2)

        plt.clf()


plot_data()
