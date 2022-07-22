import sys
import numpy as np
import matplotlib.pyplot as plt

plt.style.use('../../plot_scripts/plot_style.mplstyle')

# sys.path.insert(0, '../')
import para_scan as ps

ncells = 64*64*640

col_list    = ['tab:red', 'tab:orange', 'tab:green', 'tab:blue']
legend_list = ['B_x'    , 'B_y'       , 'B_z'      , 'hydro'   ]

for i in range(ps.amb_rho):
    for j in range(ps.Ma):

        plt.figure()

        for k in range(3):
            for B_fl in [True, False]:

                if not(B_fl) and k!=0:
                    continue

                if not(B_fl):
                    plot_i = 3
                else:
                    plot_i = k
                
                file_add = ps.filename_mix_add_ext(i,j,k,B_fl)
                dir_name = f'mix{file_add}'

                print("Accessing {dirname} ...")

                try:
                    luminosity   = np.loadtxt(f"save_arr/{dir_name}luminosity{file_add}")
                    time_list    = np.loadtxt(f"save_arr/{dir_name}time_list{file_add}")
                    B            = np.loadtxt(f"save_arr/{dir_name}B{file_add}")
                    entanglement = np.loadtxt(f"save_arr/{dir_name}entanglement{file_add}")
                except:
                    continue


                if B_fl:

                    B_x = B[:,0]
                    B_y = B[:,1]
                    B_z = B[:,2]

                    B_mag = np.sqrt(B_x*B_x + B_y*B_y + B_z*B_z)

                    # plt.plot(time_list, 0.17/B_mag, color=col_list[i_fl], label=legend_list[i_fl])
                    # plt.plot(time_list, B_z, color=col_list[i_fl], label=legend_list[i_fl])
                    # plt.plot(time_list, entanglement, color=col_list[plot_i], label=legend_list[plot_i])

                plt.plot(time_list, luminosity, color=col_list[plot_i], label=legend_list[plot_i])

plt.yscale('log')
plt.legend(loc='lower right')

plt.xlabel('time (Myr)')
plt.ylabel('Luminosity')
# plt.ylim(1e3, 1e9)
plt.ylabel(r'$\delta $B$_{\rm z, rms}$')

plt.axhline(0.037, linestyle='dashed')