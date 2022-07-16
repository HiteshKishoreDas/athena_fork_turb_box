#%%

import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sg

import read_tab as rt

plt.style.use('../plot_scripts/plot_style.mplstyle')

def dA (A):
    return np.roll(A,-1) - A


coord_dir  = ["","cartesian_1d/"]
coord_list = ["spherical","cartesian"]

# coord_dir  = ["cartesian_1d/"]
# coord_list = ["cartesian"]

# coord_dir  = [""]
# coord_list = ["spherical"]




fig1, ax1 = plt.subplots(1,1,figsize=(15,10))
fig2, ax2 = plt.subplots(1,1,figsize=(15,10))

for i_c, coord in enumerate(coord_dir):

    R_in_list = []
    R_ou_list = []
    time_list = []

    for N in range(1,401):

        print(f'Reading file {N}...')

        fn = f'{coord}Sod.block0.out1.{str(N).zfill(5)}.tab'
        prim, time_sim = rt.read_tab(fn)

        x  = prim['x1v']
        dx = dA(x)

        rho = prim['rho']
        drho = dA(rho)/dx 

        R_in = x[np.argmax(drho)]
        R_ou = x[np.argmin(drho)]

        R_in_list.append(R_in)
        R_ou_list.append(R_ou)
        time_list.append(time_sim)

        if time_sim>1.6:
            break

    R_in_list = sg.savgol_filter(R_in_list, window_length=11, polyorder=3)
    R_ou_list = sg.savgol_filter(R_ou_list, window_length=11, polyorder=3)

    vel_in = dA(R_in_list)/dA(time_list)
    vel_ou = dA(R_ou_list)/dA(time_list)

    vel_in = sg.savgol_filter(vel_in, window_length=101, polyorder=3)
    vel_ou = sg.savgol_filter(vel_ou, window_length=101, polyorder=3)

    ax1.plot(time_list, R_in_list-R_in_list[0]+0.5, label=f"Inner edge: {coord_list[i_c]}")
    ax1.plot(time_list, R_ou_list-R_ou_list[0]+0.5, label=f"Outer edge: {coord_list[i_c]}")

    ax2.plot(time_list, vel_in, label=f"Inner edge: {coord_list[i_c]}")
    ax2.plot(time_list, vel_ou, label=f"Outer edge: {coord_list[i_c]}")

plt.legend()

ax1.set_xlim(0,1.6)
ax2.set_xlim(0,1.6)

ax1.set_ylim(0,0.55)
# ax2.set_ylim(-1,1)

# %%
