#%%

import numpy as np
import yt
from yt import derived_field
# import h5py as h5
from globals import g


# Temperature in Kelvin
# @derived_field(name="temp", units="code_length**3*code_pressure/code_mass", \
#     sampling_type="cell",force_override=True)
# def _temp(field, data):
#     return (data["gas", "pressure"] / data["gas", "density"]) * KELVIN * mu

#%%

## To plot slices or projections ##
#_________________________________#

def slice_plot (ds, z_range=[0,0], z_lim_flag=False, field='density', axis='x'):

    p = yt.SlicePlot(ds,axis,("gas",field))

    p.annotate_title(f"{field} slice")

    if (z_lim_flag):
        p.set_zlim(("gas", field), z_range[0], z_range[1])

    return p

    # p.save(f'Plots/Slices/{field}_{N}.png')

def projection_plot (ds,field='density', wt_field='density', axis='x'):

    p = yt.ProjectionPlot(ds,axis,("gas",field), weight_field=("gas",wt_field))

    p.annotate_title("Temperature slice")

    return p

    # p.save(f'Plots/Projections/{field}_wt{wt_field}_{N}.png')


#%%

## To make slice and projection videos ##
# _______________________#

from matplotlib import rc_context
from matplotlib.animation import FuncAnimation

def slice_vid (ts, z_range=[0,0], z_lim_flag=False, field='density', axis='x',cmap='hot'):

    plot = yt.SlicePlot(ts[0], axis, ("gas", field))

    fig = plot.plots[("gas", field)].figure

    def animate(i):
        ds = ts[i]
        plot._switch_ds(ds)
        plot.annotate_title(f"{field}: t = {i}")
        plot.set_cmap(field,cmap=cmap)

        if (z_lim_flag):
            plot.set_zlim(("gas", field), z_range[0], z_range[1])

        # plot.annotate_velocity(factor=16)

        # print(i)

    animation = FuncAnimation(fig, animate, frames=len(ts),interval=100)

    # # Override matplotlib's defaults to get a nicer looking font
    # with rc_context({"mathtext.fontset": "stix"}):
    #     animation.save(vidpath)


    print("____________________________________________________")
    print("Video made")

    return animation

def projection_vid (ts, z_range=[0,0], z_lim_flag=False, field='density', wt_field='density', axis='x',cmap='hot'):

    plot = yt.ProjectionPlot(ts[0], axis, ("gas", field),weight_field=("gas", wt_field))

    fig = plot.plots[("gas", field)].figure

    def animate(i):
        ds = ts[i]
        plot._switch_ds(ds)
        plot.annotate_title(f"{field}: t = {i}")
        plot.set_cmap(field,cmap=cmap)

        if (z_lim_flag):
            plot.set_zlim(("gas", field), z_range[0], z_range[1])
            
        # plot.annotate_velocity(factor=16)

        # print(i)

    animation = FuncAnimation(fig, animate, frames=len(ts),interval=100)

    # # Override matplotlib's defaults to get a nicer looking font
    # with rc_context({"mathtext.fontset": "stix"}):
    #     animation.save(vidpath)


    print("____________________________________________________")
    print("Video made")

    return animation



#%%

## To plot phase plots ##
#_______________________#

# def phase_plot(ds):
#     from yt import derived_field

#     data_dir = "../work_turb/freya_output/128^3_cells_64_cores/"

#     @derived_field(name="mass", units="g", sampling_type="cell")
#     def _mass(field, data):
#         return data["gas", "density"] * data["athena_pp", "cell_volume"]

#     @derived_field(name="temp", units="code_length**3*code_pressure/code_mass", sampling_type="cell")
#     def _temp(field, data):
#         return data["gas", "pressure"] / data["gas", "density"]

#     # ds = yt.load(data_dir+"Turb.out2.00020.athdf")

#     my_sphere = ds.sphere("c", (0.4, "cm"))

#     plot = yt.PhasePlot(
#         my_sphere,
#         ("gas", "density"),
#         ("gas", "temp"),
#         ("gas", "mass"),
#         weight_field=None,
#     )

#     return plot


# data_dir = "../work_turb/freya_output/128^3_cells_64_cores/"
# ds0 = yt.load(data_dir+"Turb.out2.00020.athdf")
# phase_plot(ds0).show()



#%%

# ## To plot phase plot video ##
# #_______________________#

# from matplotlib import rc_context
# from matplotlib.animation import FuncAnimation

# data_dir = "../work_turb/freya_output/128^3_cells_64_cores/"

# ts = yt.load(data_dir+"Turb.out2.*.athdf")

# plot = phase_plot(ts[0])
# # plot = yt.SlicePlot(ts[0], "z", ("gas", "density"))
# # plot.set_zlim(("gas", "density"), 8e-29, 3e-26)

# fig = plot.plots[("gas", "mass")].figure

# # animate must accept an integer frame number. We use the frame number
# # to identify which dataset in the time series we want to load
# def animate(i):
#     ds = ts[i]
#     plot._switch_ds(ds)
#     plot.annotate_title("Phase plot: t = "+str(i))
#     # plot.annotate_velocity(factor=16)


# animation = FuncAnimation(fig, animate, frames=len(ts),interval=100)

# # Override matplotlib's defaults to get a nicer looking font
# with rc_context({"mathtext.fontset": "stix"}):
#     animation.save("phase_plot_mov.mp4")
# %%
