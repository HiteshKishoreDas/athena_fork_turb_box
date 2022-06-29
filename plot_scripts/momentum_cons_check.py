import numpy as np
import yt
import h5py as h5

from matplotlib import pyplot as plt

data_dir = "../work_turb_cloud/"

den_mid = []
vel_mid1 = []
vel_mid2 = []
vel_mid3 = []
times = []

ts = yt.load(data_dir+"Turb.out2.*.athdf")
# ds = yt.load(data_dir+"Turb.out2.00010.athdf")

for ds in ts:

    times.append(ds.current_time)

    all_data_level_0 = ds.covering_grid(
        level=0, left_edge=[0, 0.0, 0.0], dims=ds.domain_dimensions
    )

    den_mid.append(np.array(all_data_level_0["rho"])[32,32,32])
    vel_mid1.append(np.array(all_data_level_0["vel1"])[32,32,32])
    vel_mid2.append(np.array(all_data_level_0["vel2"])[32,32,32])
    vel_mid3.append(np.array(all_data_level_0["vel3"])[32,32,32])


#%%

den_mid  = np.array(den_mid)
vel_mid1 = np.array(vel_mid1)
vel_mid2 = np.array(vel_mid2)
vel_mid3 = np.array(vel_mid3)
times = np.array(times)


plt.figure(figsize=(8,8))
# plt.xlim(0.45,0.6)
# plt.plot(times,vel_mid1)
# plt.plot(times,vel_mid2)
# plt.plot(times,vel_mid3)
plt.plot(times,den_mid*(vel_mid1**2+vel_mid2**2+vel_mid3**2))
# plt.plot(times,vel_mid1)
# %%
