import matplotlib.pyplot as plt
from math import fabs, sqrt
import seaborn as sb
import numpy as np
from matplotlib import cm
from matplotlib.ticker import LinearLocator

sb.set_style('whitegrid')

# Plot the surface.
def show_3d_graph(X, Y, Z, center=None, test=None):
    np.set_printoptions(precision=0)
    X, Y = np.meshgrid(X, Y, indexing='ij')

    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

    # plot_wireframe, plot_surface
    surf = ax.plot_surface(X, Y, Z, cmap=cm.inferno, \
    edgecolor='darkred', linewidth=0.1)

    # cset = ax.contourf(X, Y, Z, zdir='z', cmap=cm.inferno)

    # Customize the z axis.
    #ax.set_zlim(-1.01, 1.01)
    ax.zaxis.set_major_locator(LinearLocator(10))
    # A StrMethodFormatter is used automatically
    ax.zaxis.set_major_formatter('{x:.02f}')

    # Add a color bar which maps values to colors.
    fig.colorbar(surf, shrink=0.5, aspect=5)

    if center:
        x, y = center
        z = np.nanmax(Z) * 1.01
        ax.scatter(x, y, z, marker='o')

    if test:
        x, y, z = test
        ax.scatter(x, y, z, marker='^')

    #plt.axis('equal')
    #ax.set_box_aspect(aspect = (1, 1, 1))
    plt.show()
