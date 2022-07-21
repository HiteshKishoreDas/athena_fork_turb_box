import sys
import numpy as np
import matplotlib.pyplot as plt

plt.style.use('../plot_scripts/plot_style.mplstyle')


dir_name = "low_beta/"

file_list = ['test_neun_hyd_cool','test_sechs_Bx_cool','test_seben_By_cool','test_acht_Bz_cool']
# file_list = ['test_sechs_Bx_cool']#,'test_seben_By_cool','test_acht_Bz_cool']

ncells = 64*64*640

col_list = ['tab:blue', 'tab:red', 'tab:orange', 'tab:green']
legend_list = ['hydro', 'B_x', 'B_y', 'B_z']

plt.figure()

for i_fl,fl in enumerate(file_list):

    print(f'i_fl: {i_fl}')

    luminosity   = np.loadtxt(f"save_arr/{dir_name}luminosity_{fl}")
    time_list    = np.loadtxt(f"save_arr/{dir_name}time_list_{fl}")
    B            = np.loadtxt(f"save_arr/{dir_name}B_{fl}")
    entanglement = np.loadtxt(f"save_arr/{dir_name}entanglement_{fl}")

    if i_fl!=0:
        B_x = B[:,0]
        B_y = B[:,1]
        B_z = B[:,2]

        B_mag = np.sqrt(B_x*B_x + B_y*B_y + B_z*B_z)

        # plt.plot(time_list, 0.17/B_mag, color=col_list[i_fl], label=legend_list[i_fl])
        # plt.plot(time_list, B_z, color=col_list[i_fl], label=legend_list[i_fl])
        plt.plot(time_list, entanglement, color=col_list[i_fl], label=legend_list[i_fl])

    # plt.plot(time_list, luminosity, color=col_list[i_fl], label=legend_list[i_fl])

plt.yscale('log')
plt.legend(loc='lower right')

plt.xlabel('time (Myr)')
plt.ylabel('Luminosity')
# plt.ylim(1e3, 1e9)
plt.ylabel(r'$\delta $B$_{\rm z, rms}$')

plt.axhline(0.037, linestyle='dashed')