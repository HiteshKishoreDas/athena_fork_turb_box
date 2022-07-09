import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import yt

# plt.style.use('../plot_scripts/plot_style.mplstyle')

import para_scan as ps
from globals import *


fig = plt.figure()

N_list = [501,502,530]

def filename(N):
    file_name = f"para_scan{ps.filename_cloud_add(2,0)}/Turb.out2.00{N}.athdf"

    return file_name

fn1a = "para_scan_Rlsh0_100_res0_128_rseed_1_M_0.5_chi_100_hydro/Turb.out2.00501.athdf"
fn2a = "para_scan_Rlsh0_100_res0_128_rseed_1_M_0.5_chi_100_hydro/Turb.out2.00510.athdf"
fn3a = "para_scan_Rlsh0_100_res0_128_rseed_1_M_0.5_chi_100_hydro/Turb.out2.00565.athdf"

fn1b = "para_scan_Rlsh0_100_res0_128_rseed_1_M_0.5_chi_100_beta_100/Turb.out2.00501.athdf"
fn2b = "para_scan_Rlsh0_100_res0_128_rseed_1_M_0.5_chi_100_beta_100/Turb.out2.00510.athdf"
fn3b = "para_scan_Rlsh0_100_res0_128_rseed_1_M_0.5_chi_100_beta_100/Turb.out2.00565.athdf"


fns = [fn1a, fn2a, fn3a, fn1b, fn2b, fn3b]

grid = AxesGrid(
    fig,
    (0.075, 0.075, 0.85, 0.85),
    nrows_ncols=(2, 3),
    axes_pad=0.05,
    label_mode="L",
    share_all=True,
    cbar_location="right",
    cbar_mode="single",
    cbar_size="3%",
    cbar_pad="0%",
)

units_override = {
    "length_unit": (1.0, "kpc"),
    "time_unit": (1.0, "Myr"),
    "mass_unit": (2.454e7, "Msun"),
}

for i, fn in enumerate(fns):
    # Load the data and create a single plot
    ds = yt.load(fn, units_override=units_override)  # load data
    p = yt.ProjectionPlot(ds, "z", ("gas", "density"))#, weight_field="density")

    # Ensure the colorbar limits match for all plots
    # p.set_zlim(("gas", "density"), 1e0, 1e1)
    p.set_zlim(("gas", "density"), 2e-3, 1e-2)

    # This forces the ProjectionPlot to redraw itself on the AxesGrid axes.
    plot = p.plots[("gas", "density")]
    plot.figure = fig
    plot.axes = grid[i].axes
    plot.cax = grid.cbar_axes[i]

    # Finally, this actually redraws the plot.
    p._setup_plots()

p.show()
# plt.savefig("multiplot_2x2_time_series.png")