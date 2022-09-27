import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import yt

# plt.style.use('../plot_scripts/plot_style.mplstyle')

import para_scan as ps
from globals import *

units_override = {
    "length_unit": (1.0, "kpc"),
    "time_unit": (1.0, "Myr"),
    "mass_unit": (2.454e7, "Msun"),
}

for N in range(501, 601):

    fn1a = f"para_scan_Rlsh1_50_res0_256_rseed_1_M_0.5_chi_100_hydro/Turb.out2.{str(N).zfill(5)}.athdf"
    fn1b = f"para_scan_Rlsh1_50_res0_256_rseed_1_M_0.5_chi_100_beta_100/Turb.out2.{str(N).zfill(5)}.athdf"

    fig=plt.figure()

    fns = [fn1a, fn1b]

    grid = AxesGrid(
        fig,
        (0.075, 0.075, 0.85, 0.85),
        nrows_ncols=( 1, 2),
        axes_pad=0.05,
        label_mode="L",
        share_all=True,
        cbar_location="right",
        cbar_mode="single",
        cbar_size="3%",
        cbar_pad="0%",
    )


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