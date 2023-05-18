import numpy as np
from matplotlib import markers, pyplot as plt

from scaling_data import scale_neg_rseed, scale_pos_rseed

from matplotlib.patches import Patch
from matplotlib.lines import Line2D

legend_elements = [
    Line2D(
        [0],
        [0],
        marker="o",
        color="w",
        label="64^3 box, -ve rseed",
        markerfacecolor="tab:blue",
        markeredgecolor="k",
        markersize=10,
    ),
    Line2D(
        [0],
        [0],
        marker="o",
        color="w",
        label="128^3 box, -ve rseed",
        markerfacecolor="tab:red",
        markeredgecolor="k",
        markersize=10,
    ),
    Line2D(
        [0],
        [0],
        marker="D",
        color="w",
        label="64^3 box, +ve rseed",
        markerfacecolor="tab:blue",
        markeredgecolor="k",
        markersize=10,
    ),
]

# plt.show()


plt.figure()
# plt.yscale("log")
# plt.title("Strong scaling for negative rseed")
for scl_dt_neg in scale_neg_rseed:
    if scl_dt_neg.ncells == 64**3:
        plt.scatter(
            scl_dt_neg.nproc,
            scl_dt_neg.time * scl_dt_neg.nproc / scl_dt_neg.ncells * 1000,
            color="tab:blue",
            edgecolors="k",
            s=350,
        )

    elif scl_dt_neg.ncells == 128**3:
        plt.scatter(
            scl_dt_neg.nproc,
            scl_dt_neg.time * scl_dt_neg.nproc / scl_dt_neg.ncells * 1000,
            color="tab:red",
            edgecolors="k",
            s=150,
        )

    else:
        plt.scatter(
            scl_dt_neg.nproc,
            scl_dt_neg.time * scl_dt_neg.nproc / scl_dt_neg.ncells * 1000,
            color="tab:orange",
            edgecolors="k",
            s=200,
        )

for scl_dt_neg in scale_pos_rseed:
    if scl_dt_neg.ncells == 64**3:
        plt.scatter(
            scl_dt_neg.nproc,
            scl_dt_neg.time * scl_dt_neg.nproc / scl_dt_neg.ncells * 1000,
            color="tab:blue",
            edgecolors="k",
            linewidths=3,
            s=100,
            marker="D",
        )

    elif scl_dt_neg.ncells == 128**3:
        plt.scatter(
            scl_dt_neg.nproc,
            scl_dt_neg.time * scl_dt_neg.nproc / scl_dt_neg.ncells * 1000,
            color="tab:red",
            edgecolors="k",
            s=50,
            marker="D",
        )

    else:
        plt.scatter(
            scl_dt_neg.nproc,
            scl_dt_neg.time * scl_dt_neg.nproc / scl_dt_neg.ncells * 1000,
            color="tab:orange",
            edgecolors="k",
            s=100,
            marker="D",
        )

plt.xlabel("Number of processors")
plt.ylabel("Time taken (ms) * nprocs /#cells")
plt.legend(handles=legend_elements)
plt.savefig("weak_scaling.png")


plt.figure()
# plt.yscale("log")
# plt.title("Strong scaling for negative rseed")
for scl_dt_neg in scale_neg_rseed:
    if scl_dt_neg.ncells == 64**3:
        plt.scatter(
            scl_dt_neg.nproc,
            scl_dt_neg.time / 3600,
            color="tab:blue",
            edgecolors="k",
            s=300,
        )

    elif scl_dt_neg.ncells == 128**3:
        plt.scatter(
            scl_dt_neg.nproc,
            scl_dt_neg.time / 3600,
            color="tab:red",
            edgecolors="k",
            s=100,
        )

    else:
        plt.scatter(
            scl_dt_neg.nproc,
            scl_dt_neg.time / 3600,
            color="tab:orange",
            edgecolors="k",
            s=200,
        )

for scl_dt_neg in scale_pos_rseed:
    if scl_dt_neg.ncells == 64**3:
        plt.scatter(
            scl_dt_neg.nproc,
            scl_dt_neg.time / 3600,
            color="tab:blue",
            edgecolors="k",
            s=100,
            linewidths=3,
            marker="D",
        )

    elif scl_dt_neg.ncells == 128**3:
        plt.scatter(
            scl_dt_neg.nproc,
            scl_dt_neg.time / 3600,
            color="tab:red",
            edgecolors="k",
            s=50,
            marker="D",
        )

    else:
        plt.scatter(
            scl_dt_neg.nproc,
            scl_dt_neg.time / 3600,
            color="tab:orange",
            edgecolors="k",
            s=100,
            marker="D",
        )

plt.xlabel("Number of processors")
plt.ylabel("Time taken (hr)")

plt.legend(handles=legend_elements, loc="lower right")
plt.savefig("time_taken_vs_proc.png")
