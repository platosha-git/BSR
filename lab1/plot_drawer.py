import matplotlib
matplotlib.use('tkagg')

import matplotlib.pyplot as plt
from math import fabs, sqrt
import seaborn as sb
import numpy as np
from matplotlib import cm
from matplotlib.ticker import LinearLocator


def show_3d_graph(X, Y, Z, center=None, test=None):
    np.set_printoptions(precision=0)
    X, Y = np.meshgrid(X, Y, indexing='ij')

    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

    surf = ax.plot_surface(X, Y, Z, cmap=cm.inferno, edgecolor='darkred', linewidth=0.1)
    fig.colorbar(surf, shrink=0.5, aspect=5)

    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter('{x:.01f}')

    if center:
        x, y = center
        z = np.nanmax(Z) * 1.01
        ax.scatter(x, y, z, marker='o')

    plt.show()
