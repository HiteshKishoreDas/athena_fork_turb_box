import sys
import numpy as np
import matplotlib.pyplot as plt

import scipy.signal as sg

# plt.style.use('../plot_scripts/plot_style.mplstyle')
plt.style.use('dark_background')

# sys.path.insert(0, '../')
sys.path.insert(0, './plot_scripts')
import para_scan as ps
import hst_read as hr

import globals as g

ncells = 64*64*640

col_list    = ['tab:red', 'tab:orange', 'tab:green', 'tab:blue']
legend_list = ['B_x'    , 'B_y'       , 'B_z'      , 'hydro'   ]


def lum_fn(hst,i):

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
    dz = (ps.box_length[i]/ps.nx3[0])
    dy = (ps.box_width[i] /ps.nx2[0])
    dx = (ps.box_width[i] /ps.nx1[0])

    luminosity  = (dcool/dt)[1:-1] * dx * dy * dz
    luminosity /= (ps.box_width[i] * ps.box_width[i])
    time       = time[1:-1]

    return time, luminosity
    # return hst.time, tot_cool 


# fig, ax = plt.subplots(nrows=2, ncols=2, figsize = (20,25))
# fig, ax = plt.subplots(nrows=len(ps.box_width), ncols=1, figsize = (10,85))
fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (10,10))

Da_list       = []
L_avg_list    = []
Q_theory_list = []

label_list=['Large box width', r'Large $\Lambda_0$']

Lambda_fac = [ 1.0, 1.0, 1.0, 1.0, 1.0, 10.0, 100.0, 1000.0, 10000.0]

# file_list = ['cooling_test', 'cooling_test_2']
# file_list = ['hst_test']

for i_plot,i in enumerate(range(len(ps.box_width))):
    for j in range(len(ps.Ma)):
        for k in range(1):
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
                file_add = ps.filename_mix_add_ext(i,j,k,B_fl)# [:-7]
                if not B_fl:
                    setup_name = file_add
                dir_name = f'mix{file_add}'

                # dir_name = file_list[i] 

                try:
                    hst = hr.hst_data(f'{dir_name}/Turb.hst', ncells, B_fl)     
                except:
                    continue

                time, luminosity = lum_fn(hst, i)

                # ax[i_plot].plot(time/ps.t_KH[i_plot], luminosity)#, color=col_list[plot_i], label=legend_list[plot_i])
                # ax.plot(time/ps.t_KH[i_plot], luminosity, label=label_list[i_plot])#, color=col_list[plot_i], label=legend_list[plot_i])
                ax.plot(time/ps.t_KH[i_plot], luminosity, label=f'Lbox = {ps.box_width[i_plot]}, '+ r'$\Lambda_0$=' + f'{Lambda_fac[i_plot]}')#, color=col_list[plot_i], label=legend_list[plot_i])
                # ax.plot(time, luminosity)#, color=col_list[plot_i], label=legend_list[plot_i])





                if not B_fl:

                    # L_avg = np.average(luminosity[-1250:-1000])
                    L_avg = np.average(luminosity[-250:])
                    print(f'L_avg: {L_avg}')

                    # ax[i_plot].axhline(L_avg, linestyle='dashed')#,\
                    # ax.axhline(L_avg, linestyle='dashed', label=f'L_avg = {"%.4e" % L_avg}')
                          # color=col_list[plot_i],\
                        # label=r'L$_{\rm avg}$'+ f' = {np.round(L_avg,3)}')


        # ax[i].set_yscale('log')
        # ax[i].legend(loc='lower right')

        # ax[i].set_xlabel('time (Myr)')
        # ax[i].set_ylabel('Luminosity')

        # ax[i_plot].set_yscale('log')
        ax.set_yscale('log')
        # ax[i].legend(loc='lower right')

        # ax[i].set_title(f'Lambda_frac: {ps.Lambda_fac[i]}')

        t_turb = ps.box_width[i_plot]/ps.v_shear
        t_cool = ps.t_cool_Da

        if i_plot==0: 
            Da = (t_turb/t_cool)[0]

        # ax[i_plot].set_title(f'Box_width: {ps.box_width[i_plot]} kpc, Da = {Da}')
        # ax[i_plot].set_xlabel('time (Myr)')
        # ax[i_plot].set_ylabel('Luminosity')

        print(f'Box_width: {ps.box_width[i_plot]} kpc, Da = {Da}')

        # ax.set_title(f'Box_width: {ps.box_width[0]} kpc, Da = {Da}')
        # ax.set_title(f'Da = {np.round(Da,3)}')
        ax.set_xlabel('time (Myr)')
        ax.set_ylabel('Luminosity')
        # ax[i,j].set_ylim(0.99,1)

        # ax[i].set_title(r"$\mathcal{M}_{\rm a}$"+f" = {ps.Ma[j]}; " + r"t$_{\rm cool}$/t$_{\rm KH}$" + f" = {tcool_tKH}" + f"; Da = {np.round(Da, 3)}")#+f"\n {setup_name}")
        
        plt.legend(fontsize=12)
# plt.savefig("test.png")

#%%
