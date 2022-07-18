#%%

import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sg

import read_tab as rt
import mass_cons_check as mc

plt.style.use('../plot_scripts/plot_style.mplstyle')

def dA (A):
    return np.roll(A,-1) - A


coord_dir  = ["","cartesian_1d/", "larger_R/"]
coord_list = ["spherical","cartesian", "larger_R"]

R_list = [0.5,0.5,5]
t_cut_list = [1.6,1.6,16]

# coord_dir  = ["cartesian_1d/"]
# coord_list = ["cartesian"]

# coord_dir  = [""]
# coord_list = ["spherical"]

color_list = ['tab:orange', 'tab:green', 'tab:blue']
lines_list = ['solid'   , 'dashed' ]
line_size  = [8,5,3]

fig1, ax1 = plt.subplots(1,1,figsize=(15,10))
fig2, ax2 = plt.subplots(1,1,figsize=(15,10))

for i_c, coord in enumerate(coord_dir):

    R_in_list = []
    R_ou_list = []
    time_list = []

    print(i_c)

    for N in range(1,401):

        # print(f'Reading file {N}...')

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

        # if time_sim>t_cut_list[i_c]:
        #     break

    R_in_list = np.array(R_in_list)
    R_ou_list = np.array(R_ou_list)
    time_list = np.array(time_list)

    R_in_list = sg.savgol_filter(R_in_list, window_length=11, polyorder=3)
    R_ou_list = sg.savgol_filter(R_ou_list, window_length=11, polyorder=3)

    vel_in = dA(R_in_list)/dA(time_list)
    vel_ou = dA(R_ou_list)/dA(time_list)

    vel_in = sg.savgol_filter(vel_in, window_length=101, polyorder=3)
    vel_ou = sg.savgol_filter(vel_ou, window_length=101, polyorder=3)

    y_plot_in = (R_in_list-R_in_list[0]+R_list[i_c])/R_list[i_c]
    y_plot_ou = (R_ou_list-R_ou_list[0]+R_list[i_c])/R_list[i_c]

    ax1.plot(time_list/R_list[i_c], y_plot_in, label=f"Inner edge: {coord_list[i_c]}", color=color_list[i_c], linestyle='dashed', linewidth=line_size[i_c])
    ax1.plot(time_list/R_list[i_c], y_plot_ou, label=f"Outer edge: {coord_list[i_c]}", color=color_list[i_c], linestyle='solid' , linewidth=line_size[i_c])
    
    # ax1.plot(time_list, R_list[i_c] - mc.u_s*time_list, label="Incident shock theory", linestyle="dotted", color='tab:red' )
    # ax1.plot(mc.t_arr , R_list[i_c] - mc.R_Mcons      , label="Outer shell theory"   , linestyle="dotted", color='tab:blue')

    ax2.plot(time_list/R_list[i_c], vel_in, label=f"Inner edge: {coord_list[i_c]}", color=color_list[i_c], linestyle='dashed', linewidth=line_size[i_c])
    ax2.plot(time_list/R_list[i_c], vel_ou, label=f"Outer edge: {coord_list[i_c]}", color=color_list[i_c], linestyle='solid' , linewidth=line_size[i_c])

    # ax2.axhline(-mc.u_s, label=f"Incident shock theory", color='tab:red', linestyle="dotted")

ax1.legend(loc='upper left')
ax2.legend(loc='upper left')

ax2.axhline(0.0, color='k', linestyle='dotted')

# ax2.axhline(mc.u_s, color='k', linestyle='-.')
# ax2.axhline(-mc.u_s, color='k', linestyle='-.')

ax1.set_xlabel('t/R')
ax1.set_ylabel(r'R/R$_0$')

ax2.set_xlabel('t/R')
ax2.set_ylabel('v')

# ax1.set_xlim(0,1.6)
# ax2.set_xlim(0,1.6)

# ax1.set_ylim(0,1.1)
ax2.set_ylim(-0.75,0.5)

# %%
