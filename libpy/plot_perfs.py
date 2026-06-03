import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("perfs",skiprows=1)

x      = data[:, 0]   # col 1: bar positions
widths = data[:, 1]   # col 2: bar widths
cols   = [2, 3, 4, 5, 6, 8]  # 0-indexed cols 3,4,5,6,7,9
labels = ["mgp.i", "mpg.nrelax", "n_grid", "perf.t (s)", "perf.speed", "perfs.ispeed"]

fig, axes = plt.subplots(2, 3, figsize=(14, 7))

for ax, col, label in zip(axes.flat, cols, labels):
    ax.bar(x, data[:, col], width=widths, align="edge")
    ax.set_xlabel("time")
    ax.set_ylabel(label)

plt.tight_layout()
fig.savefig('perfs.png')
plt.show()
