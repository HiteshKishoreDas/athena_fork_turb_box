import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import yt

# plt.style.use('../plot_scripts/plot_style.mplstyle')

import para_scan as ps
from globals import *
import v_turb as vt


M = 0.25
v_turb_predict = M*vt.cs_calc(ps.T_hot_req,ps.mu)
t_eddy = ps.L_box/v_turb_predict

fig = plt.figure()

N_list = [501,510,565]

R_lsh = ["5","50"]

def filename(N):
    file_name = f"para_scan{ps.filename_cloud_add(2,0)}/Turb.out2.00{N}.athdf"

    return file_name

fn1a = f"para_scan_Rlsh{R_lsh[0]}_{R_lsh[1]}_res0_128_rseed_1_M_{M}_chi_100_hydro/Turb.out2.{str(N_list[0]).zfill(5)}.athdf"
fn2a = f"para_scan_Rlsh{R_lsh[0]}_{R_lsh[1]}_res0_128_rseed_1_M_{M}_chi_100_hydro/Turb.out2.{str(N_list[1]).zfill(5)}.athdf"
fn3a = f"para_scan_Rlsh{R_lsh[0]}_{R_lsh[1]}_res0_128_rseed_1_M_{M}_chi_100_hydro/Turb.out2.{str(N_list[2]).zfill(5)}.athdf"

fn1b = f"para_scan_Rlsh{R_lsh[0]}_{R_lsh[1]}_res0_128_rseed_1_M_{M}_chi_100_beta_100/Turb.out2.{str(N_list[0]).zfill(5)}.athdf"
fn2b = f"para_scan_Rlsh{R_lsh[0]}_{R_lsh[1]}_res0_128_rseed_1_M_{M}_chi_100_beta_100/Turb.out2.{str(N_list[1]).zfill(5)}.athdf"
fn3b = f"para_scan_Rlsh{R_lsh[0]}_{R_lsh[1]}_res0_128_rseed_1_M_{M}_chi_100_beta_100/Turb.out2.{str(N_list[2]).zfill(5)}.athdf"


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

current_time = []

for i, fn in enumerate(fns):
    # Load the data and create a single plot
    ds = yt.load(fn)#, units_override=units_override)  # load data
    p = yt.ProjectionPlot(ds, "z", ("gas", "density"))#, weight_field="density")

    current_time.append(float(ds.current_time))

    # Ensure the colorbar limits match for all plots
    # p.set_zlim(("gas", "density"), 5e-2, 5e-1)
    # p.set_zlim(("gas", "density"), 2e-3, 1e-2)

    # This forces the ProjectionPlot to redraw itself on the AxesGrid axes.
    plot = p.plots[("gas", "density")]
    plot.figure = fig
    plot.axes = grid[i].axes
    plot.cax = grid.cbar_axes[i]

    # Finally, this actually redraws the plot.
    p._setup_plots()

p.show()
# plt.savefig("multiplot_2x2_time_series.png")