import sys
import numpy as np
import matplotlib.pyplot as plt

import scipy.signal as sg

plt.style.use('../plot_scripts/plot_style.mplstyle')

# sys.path.insert(0, '../')
sys.path.insert(0, './plot_scripts')
import para_scan as ps
import hst_read as hr

import globals as g

ncells = 64*64*640

col_list    = ['tab:red', 'tab:orange', 'tab:green', 'tab:blue']
legend_list = ['B_x'    , 'B_y'       , 'B_z'      , 'hydro'   ]


def lum_fn(hst):

    N_last = 2000

    #* Copied cumulative cooling data
    tot_cool = np.copy(hst.total_cooling)[:N_last]

    dcool = np.roll(tot_cool, -1) - tot_cool
    dt    = np.roll(hst.time[:N_last], -1) - hst.time[:N_last]
    # time  = hst.time

    dcool = dcool   [hst.cold_gas_fraction[:N_last]>0.1]
    dt    = dt      [hst.cold_gas_fraction[:N_last]>0.1]
    time  = (hst.time[:N_last])[hst.cold_gas_fraction[:N_last]>0.1]

    box_full = np.argwhere(hst.cold_gas_fraction[:N_last]>0.998)

    if len(box_full)!=0:
        dcool = dcool   [:np.min(box_full)]
        dt    = dt      [:np.min(box_full)]
        time  = hst.time[:np.min(box_full)]

    box_volume = ps.box_length*ps.box_width*ps.box_width
    dz = (ps.box_length[0]/ps.nx3[0])
    dy = (ps.box_width[0] /ps.nx2[0])
    dx = (ps.box_width[0] /ps.nx1[0])

    luminosity  = (dcool/dt)[1:-1] * dx * dy * dz
    luminosity /= (ps.box_width[0] * ps.box_width[0])
    time       = time[1:-1]

    return time, luminosity
    # return hst.time, tot_cool 


# fig, ax = plt.subplots(nrows=2, ncols=2, figsize = (20,25))
fig, ax = plt.subplots(nrows=len(ps.Lambda_fac), ncols=1, figsize = (10,85))
# fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (10,15))

Da_list       = []
L_avg_list    = []
Q_theory_list = []

# file_list = ['cooling_test', 'cooling_test_2']
# file_list = ['hst_test']

for i in range(len(ps.Lambda_fac)):
    for j in range(len(ps.Ma)):
        for k in range(3):
            # for B_fl in [True, False]:
            for B_fl in [False]:

                # #* Non-magnetic cases
                # if not(B_fl) and k!=0:
                #     continue

                # # For colors and linestyles
                # if not(B_fl):
                #     plot_i = 3
                # else:
                #     plot_i = k


                #* Read history file 
                file_add = ps.filename_mix_add_ext(i,j,k,B_fl)
                if not B_fl:
                    setup_name = file_add
                dir_name = f'mix{file_add}'

                # dir_name = file_list[i] 

                hst = hr.hst_data(f'{dir_name}/Turb.hst', ncells, B_fl)     
                time, luminosity = lum_fn(hst)

                ax[i].plot(time, luminosity)#, color=col_list[plot_i], label=legend_list[plot_i])
                # ax.plot(time, luminosity)#, color=col_list[plot_i], label=legend_list[plot_i])





                if not B_fl:
                    
                    L_avg = np.average(luminosity[-150:])
                    ax[i].axhline(L_avg, linestyle='dashed')#,\
                          # color=col_list[plot_i],\
                        # label=r'L$_{\rm avg}$'+ f' = {np.round(L_avg,3)}')


        # ax[i].set_yscale('log')
        # ax[i].legend(loc='lower right')

        # ax[i].set_xlabel('time (Myr)')
        # ax[i].set_ylabel('Luminosity')

        ax[i].set_yscale('log')
        # ax[i].legend(loc='lower right')

        ax[i].set_title(f'Lambda_frac: {ps.Lambda_fac[i]}')
        ax[i].set_xlabel('time (Myr)')
        ax[i].set_ylabel('Luminosity')

        # ax[i,j].set_ylim(0.99,1)

        # ax[i].set_title(r"$\mathcal{M}_{\rm a}$"+f" = {ps.Ma[j]}; " + r"t$_{\rm cool}$/t$_{\rm KH}$" + f" = {tcool_tKH}" + f"; Da = {np.round(Da, 3)}")#+f"\n {setup_name}")


#%%
