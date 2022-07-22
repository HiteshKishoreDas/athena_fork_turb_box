import sys
import numpy as np
import matplotlib.pyplot as plt

import scipy.signal as sg

plt.style.use('../plot_scripts/plot_style.mplstyle')

# sys.path.insert(0, '../')
sys.path.insert(0, './plot_scripts')
import para_scan as ps
import hst_read as hr

ncells = 64*64*640

col_list    = ['tab:red', 'tab:orange', 'tab:green', 'tab:blue']
legend_list = ['B_x'    , 'B_y'       , 'B_z'      , 'hydro'   ]

fig, ax = plt.subplots(nrows=2, ncols=2, figsize = (20,25))

for i in range(len(ps.amb_rho)):
    for j in range(len(ps.Ma)):


        for k in range(3):
            for B_fl in [True, False]:
            # for B_fl in [False]:

                if not(B_fl) and k!=0:
                    continue

                if not(B_fl):
                    plot_i = 3
                else:
                    plot_i = k
                
                file_add = ps.filename_mix_add_ext(i,j,k,B_fl)
                dir_name = f'mix{file_add}'

                hst = hr.hst_data(f'{dir_name}/Turb.hst', ncells, B_fl)     

                # try:
                #     luminosity   = np.loadtxt(f"save_arr/luminosity{file_add}")
                #     time_list    = np.loadtxt(f"save_arr/time_list{file_add}")
                #     B            = np.loadtxt(f"save_arr/B{file_add}")
                #     entanglement = np.loadtxt(f"save_arr/entanglement{file_add}")
                #     print(f"Accessing {dir_name} ...")
                # except:
                #     continue

                # ax[i,j].plot(time_list, np.abs(luminosity), color=col_list[plot_i], label=legend_list[plot_i])

                # if B_fl:

                #     B_x = B[:,0]
                #     B_y = B[:,1]
                #     B_z = B[:,2]

                #     B_mag = np.sqrt(B_x*B_x + B_y*B_y + B_z*B_z)

                    # plt.plot(time_list, 0.17/B_mag, color=col_list[i_fl], label=legend_list[i_fl])
                    # plt.plot(time_list, B_z, color=col_list[i_fl], label=legend_list[i_fl])
                    # plt.plot(time_list, entanglement, color=col_list[plot_i], label=legend_list[plot_i])


                tot_cool = np.copy(hst.total_cooling)
                # tot_cool = sg.savgol_filter(tot_cool, window_length=21, polyorder=3)

                dcool = np.roll(tot_cool, -1) - tot_cool
                dt    = np.roll(hst.time, -1) - hst.time

                luminosity = (dcool/dt)[1:-1]/(ps.box_width**2)
                time       = np.copy(hst.time[1:-1])

                ax[i,j].plot(time, luminosity, color=col_list[plot_i], label=legend_list[plot_i])

                # if B_fl:
                #     ax[i]
                # ax[i,j].plot(hst.time, np.abs(hst.total_cooling), color=col_list[plot_i], label=legend_list[plot_i])

                if not B_fl:
                    L_avg = np.average(luminosity[-100:])
                    ax[i,j].axhline(L_avg, linestyle='dashed', color=col_list[plot_i],\
                        label=r'L$_{\rm avg}$'+ f' = {np.round(L_avg,3)}')


        ax[i,j].set_yscale('log')
        ax[i,j].legend(loc='lower right')

        ax[i,j].set_xlabel('time (Myr)')
        ax[i,j].set_ylabel('Luminosity')

        tcool_tKH = np.round(ps.t_cool_mix[i]/ps.t_KH[0], 3)
        ax[i,j].set_title(r"$\mathcal{M}_{\rm a}$"+f" = {ps.Ma[j]}; " + r"t$_{\rm cool}$/t$_{\rm KH}$" + f" = {tcool_tKH}" + f"; Da = {np.round(1/(np.sqrt(ps.chi_cold)*tcool_tKH), 3)}")

        # if i==0:
        #     ax[i,j].set_ylim(1e4, 2e6)
        # elif i==1:
        #     ax[i,j].set_ylim(1e2, 1e5)
        plt.ylabel(r'$\delta $B$_{\rm z, rms}$')

        # plt.axhline(0.037, linestyle='dashed')

plt.show()