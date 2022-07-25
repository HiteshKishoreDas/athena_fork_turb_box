import numpy as np
import yt
from yt import derived_field

import sys
# sys.path.insert(0, '../vis/python')
# import athena_read as ar
import globals as g
import para_scan as ps


@derived_field(name="temp", sampling_type="cell")
def _temp(field, data):
    return (
        data["athena_pp", "press"]/data["gas", "density"] * g.KELVIN * g.mu
    )


@derived_field(name="total_pressure", sampling_type="cell", units="code_mass/(code_length*code_time**2)")
def _total_pressure(field, data):
    return (
        data["athena_pp", "press"]+data["gas", "magnetic_pressure"]
    )
# units_override = {"length_unit": (1.0, "kpc"),
#                   "time_unit"  : (1.0, "Myr"),
#                   "density_unit"  : (2.454e7, "Msun/kpc**3")}


N = 100

i = 1
j = 1
k = 2
B_fl = True

file_add = ps.filename_mix_add_ext(i, j, k, B_fl)
dir_name = f'mix{file_add}'
file_name = f'{dir_name}/Turb.out2.{str(N).zfill(5)}.athdf'

ds = yt.load(file_name)  # [:,:,:320]

field_list = ["rho", "vel1", "vel3"]
log_scale  = [True , False , False ]
anno_scale = [False, False , False ]

z_min  = [ ps.amb_rho[i]     , -0.1 , -0.001  ]
z_max  = [ 100*ps.amb_rho[i] ,  0.2 ,  0.002  ]

for i_fld, fld in enumerate(field_list):

    p = yt.SlicePlot(ds, axis='y', fields=fld)
    # p = yt.ProjectionPlot(ds,axis='y',fields="total_pressure")

    z_center = ps.box_length[0]*-0.05

    # p = yt.SlicePlot(ds,axis='y',fields=fld, \
    #                  center=[ps.box_width/2,ps.box_width/2,z_center], \
    #                  width=(0.5*1e-2,4.3e-3), \
    #                  origin='native')


    # p = yt.ProjectionPlot(ds,axis='y',fields="temp", \
    #                  center=[ps.box_width/2,ps.box_width/2,z_center], \
    #                  width=(1e-2,4.3e-3), \
    #                  origin='native')

    p.set_log(fld, log_scale[i_fld])

    if anno_scale[i_fld]:
        p.annotate_streamlines(("vel3"), ("vel1"))
        # p.annotate_streamlines(("Bcc3"), ("Bcc1"))

    # p.annotate_timestamp(corner="upper_left", draw_inset_box=True)

    # p.set_cmap(field=fld, cmap="CMRmap")
    p.set_cmap(field=fld, cmap="RdGy")

    # p.zoom(2)


    p.set_zlim(fld, z_min[i_fld], z_max[i_fld])


    p.save(f"Plots/{fld}_{str(N).zfill(5)}.png")

    p.show()
