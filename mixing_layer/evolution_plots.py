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
        data["athena_pp", "press"]/data["gas","density"]  * g.KELVIN * g.mu
    )

# units_override = {"length_unit": (1.0, "kpc"),
#                   "time_unit"  : (1.0, "Myr"),
#                   "density_unit"  : (2.454e7, "Msun/kpc**3")}

N = 5  

i = 0 
j = 0 
k = 2 
B_fl = True 

file_add = ps.filename_mix_add_ext(i,j,k,B_fl)
dir_name = f'mix{file_add}'
file_name = f'{dir_name}/Turb.out2.{str(N).zfill(5)}.athdf'

ds = yt.load(file_name) # [:,:,:320]


# p = yt.SlicePlot(ds,axis='y',fields="press")
# p = yt.ProjectionPlot(ds,axis='y',fields="temp")

z_center = ps.box_length[0]*-0.1

p = yt.SlicePlot(ds,axis='y',fields="temp", \
                 center=[ps.box_width/2,ps.box_width/2,z_center], \
                 width=(1e-2,4.3e-3), \
                 origin='native')
    
# p = yt.ProjectionPlot(ds,axis='y',fields="temp", \
#                  center=[ps.box_width/2,ps.box_width/2,z_center], \
#                  width=(1e-2,4.3e-3), \
#                  origin='native')

# p = yt.SlicePlot(ds,axis='y',fields="beta")

# p.annotate_streamlines(("vel3"), ("vel1"))
p.annotate_streamlines(("Bcc3"), ("Bcc1"))

# p.annotate_timestamp(corner="upper_left", draw_inset_box=True)

p.set_cmap(field="temp", cmap="CMRmap")
# p.set_cmap(field="temp", cmap="RdGy_r")
# p.set_cmap(field="beta", cmap="RdGy")
# p.set_cmap(field="vel1", cmap="RdGy")

# p.zoom(2)


p.set_zlim("temp",4e4, 4e6)
# p.set_zlim("rho",5, 50)
# p.set_zlim("vel3",-0.17, 0.17)
# p.set_zlim("vel1",-0.17, 0.17)
# p.set_zlim(('gas','mach_alfven'), 1e-1,1e1 )
# p.set_zlim("beta",1e-1, 1e1 )
# p.set_zlim("Ma",1e-1, 1e2)

# p.save("Plots/beta_By.png")

p.show()




# ds = yt.load(file_name)#, units_override=units_override)

# B  = data['Bcc1'] * data['Bcc1']
# B += data['Bcc2'] * data['Bcc2']
# B += data['Bcc3'] * data['Bcc3']
# B = np.sqrt(B)

# Ma_calc = 0.17/B