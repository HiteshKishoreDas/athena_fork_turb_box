import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import yt

# plt.style.use('../plot_scripts/plot_style.mplstyle')

import para_scan as ps
import v_turb as vt
from globals import *


fig = plt.figure()

N_list = [501,502,530]

def filename(N):
    file_name = f"para_scan{ps.filename_cloud_add(2,0)}/Turb.out2.00{N}.athdf"

    return file_name

N = 501

for N in range(501,566):

    print(f"Processing file {N}...")

    M = 0.5
    v_turb_predict = M*vt.cs_calc(ps.T_hot_req,ps.mu)
    t_eddy = ps.L_box/v_turb_predict

    Rlsh_i = 0 
    t_eddy = t_eddy[Rlsh_i]

    fn_hyd = f"para_scan_Rlsh0_100_res0_128_rseed_1_M_0.5_chi_100_hydro/Turb.out2.{str(N).zfill(5)}.athdf"
    fn_mhd = f"para_scan_Rlsh0_100_res0_128_rseed_1_M_0.5_chi_100_beta_100/Turb.out2.{str(N).zfill(5)}.athdf"


    fns = [fn_hyd, fn_mhd]
    fn_title = ["HD","MHD"]


    grid = AxesGrid(
        fig,
        (0.075, 0.075, 0.85, 0.85),
        nrows_ncols=(1, 2),
        axes_pad=0.05,
        label_mode="L",
        share_all=True,
        cbar_location="right",
        cbar_mode="single",
        cbar_size="3%",
        cbar_pad="0%",
    )

    units_override_list = {
        "length_unit": (1.0, "kpc"),
        "time_unit": (1.0, "Myr"),
        "mass_unit": (2.454e7, "Msun")
    }

    for i, fn in enumerate(fns):
        # Load the data and create a single plot
        ds = yt.load(fn, units_override=units_override_list)  # load data
        p = yt.ProjectionPlot(ds, "z", ("gas", "density"))#, weight_field="density")

        p.set_xlabel("x (kpc)")
        p.set_ylabel("y (kpc)")

        p.set_unit(("gas","density"),"Msun/kpc**2")

        # Ensure the colorbar limits match for all plots
        # p.set_zlim(("gas", "density"), 1e0, 1e1)
        p.set_zlim(("gas", "density"), 1e5, 1e6)
        # p.set_zlim(("gas", "density"), 1e6, 1e7)

        p.annotate_title(f"{fn_title[i]}") #: {round(float(ds.current_time)/t_eddy - 7, 3)} "+r"t$_{\rm eddy}$")    


        # This forces the ProjectionPlot to redraw itself on the AxesGrid axes.
        plot = p.plots[("gas", "density")]
        plot.figure = fig
        plot.axes = grid[i].axes
        plot.cax = grid.cbar_axes[i]

        # p.annotate_clear()
        # p.annotate_title(f"{fn_title[i]}: {round(float(ds.current_time)/t_eddy - 7, 3)} "+r"t$_{\rm eddy}$")    

        # Finally, this actually redraws the plot.
        p._setup_plots()

    # p.show()
    p.save(f"Plots/density_destroy/density_proj_{str(N).zfill(5)}.png")
    # plt.close()

print("Done!")