import numpy as np
import yt
from yt import derived_field

import sys
sys.path.insert(0, '../vis/python')
import athena_read

@derived_field(name="beta", sampling_type="cell")
def _beta(field, data):
    return (
        data["athena_pp", "press"]/data["gas","magnetic_pressure"] 
    )

# units_override = {"length_unit": (1.0, "kpc"),
#                   "time_unit"  : (1.0, "Myr"),
#                   "density_unit"  : (2.454e7, "Msun/kpc**3")}

# dirname = 'test_sechs_Bx_cool'
dirname = 'test_seben_By_cool'
# dirname = 'test_acht_Bz_cool/low_beta'
# dirname = 'test_neun_hyd_cool'

N = 25

file_name = dirname + f'/Turb.out2.{str(N).zfill(5)}.athdf'

ds = yt.load(file_name)

# data = athena_read.athdf(file_name)
# print(np.min(data['rho']))
# print(np.max(data['rho']))

# p = yt.SlicePlot(ds,axis='y',fields="press")
# p = yt.ProjectionPlot(ds,axis='x',fields="rho")
p = yt.SlicePlot(ds,axis='y',fields="beta")
# p = yt.SlicePlot(ds,axis='y',fields="beta")

# p.set_cmap(field="rho", cmap="CMRmap")
# p.set_cmap(field="rho", cmap="RdGy")
p.set_cmap(field="beta", cmap="RdGy")
# p.set_cmap(field="vel1", cmap="RdGy")

# p.zoom(2)

# p.set_zlim("rho",5e-2, 5e-1)
# p.set_zlim("rho",5, 50)
# p.set_zlim("vel3",-0.17, 0.17)
# p.set_zlim("vel1",-0.17, 0.17)
# p.set_zlim(('gas','mach_alfven'), 1e-1,1e1 )
p.set_zlim("beta",1e-1, 1e1 )
# p.set_zlim("Ma",1e-1, 1e2)

p.save("Plots/beta_By.png")

p.show()




# ds = yt.load(file_name)#, units_override=units_override)

# B  = data['Bcc1'] * data['Bcc1']
# B += data['Bcc2'] * data['Bcc2']
# B += data['Bcc3'] * data['Bcc3']
# B = np.sqrt(B)

# Ma_calc = 0.17/B