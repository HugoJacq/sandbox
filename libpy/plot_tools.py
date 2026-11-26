import matplotlib as mpl


def nice_slope(ax: mpl.axes.Axes, slope, X, Y0, Xt, Yt, fs, c, var):
    """
    Plot a slope on ax.
    X (tuple): move along x axis
    Y0: move vertically.
    (Xt,Yt): the position of the label, in data coordinates.
    'fs': fontsize
    'c': color of line and text
    """
    x0, x1 = X
    Y1 = Y0 * (X[1] / X[0]) ** slope
    ax.loglog([X[0], X[1]], [Y0, Y1], "--", lw=1, c=c)
    ax.text(Xt, Yt, rf"${var}^{slope}$", fontsize=fs, va="center", c=c)
