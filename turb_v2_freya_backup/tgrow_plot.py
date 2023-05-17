import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mt 
import scipy.signal as sg

plt.style.use('../plot_scripts/plot_style.mplstyle')

import para_scan as ps
import v_turb as vt
import athena_read as ar
from globals import *
import tgrow as tg

MHD = [True, False]
linestyle_list = ['solid', 'dashed']

# R_lsh = np.array([5,10,50,250,500])         # For M = 0.25
R_lsh = np.array([10,50,100,250,500,1000])     # For M = 0.5
# R_lsh = np.array([250,500,1000,2500,5000])  # For M = 0.9

MHD = [True, False]
linestyle_list = ['solid', 'dashed']

cmap_name = 'viridis' #'Paired'

cmap = mt.cm.get_cmap(cmap_name)

# cb_qnt = np.scopy(ps.dedt)
cb_qnt = np.copy(np.log10(R_lsh))

# line_col = cmap(np.log10(beta_list)/np.log10(beta_list).max())
line_col = cmap(cb_qnt/cb_qnt.max())

res = 256 #128

M = 0.5

l_sh = vt.cs_calc(ps.T_floor,ps.mu)*ps.t_cool_cloud
cloud_radius_temp = R_lsh*l_sh
L_box = cloud_radius_temp*40

v_turb_predict = M*vt.cs_calc(ps.T_hot_req,ps.mu)
t_eddy = L_box/v_turb_predict

plt.figure()

for i_MHD, MHD_flag in enumerate(MHD):
    for i in [4,5]:#range(3,len(ps.R_lsh)):

        # fn_suffix = ps.filename_cloud_func(i,j=0,rseed=1,Mach=0.5,cloud_chi=100,beta=100,MHD_flag=MHD_flag)
        if MHD_flag:
            fn_suffix = f'_Rlsh{i}_{R_lsh[i]}_res0_256_rseed_1_M_{M}_chi_100_beta_100'
        else:
            fn_suffix = f'_Rlsh{i}_{R_lsh[i]}_res0_256_rseed_1_M_{M}_chi_100_hydro'

        save_arr = np.loadtxt(f"save_arr/save_arr{fn_suffix}")

        time     = save_arr[0,:]
        cold_gas = save_arr[1,:]

        x_data = time
        y_data = cold_gas/cold_gas[0]

        # y_data= sg.savgol_filter(y_data, window_length=11, polyorder=3)

        window_hw = int(51/2)
        L_data = len(y_data)

        fit_region_y = [y_data[i-window_hw:i+window_hw] for i in range(window_hw,L_data-window_hw)]
        fit_region_x = [x_data[i-window_hw:i+window_hw] for i in range(window_hw,L_data-window_hw)]
        fit_x     = np.array([x_data[i] for i in range(window_hw,L_data-window_hw)])
        fit_y     = np.array([y_data[i] for i in range(window_hw,L_data-window_hw)])

        fit_list = [ np.polyfit(fit_region_x[i], fit_region_y[i], 1) for i in range(len(fit_x))]

        # dy = np.roll(y_data,-1) - y_data
        # dx = np.roll(x_data,-1) - x_data

        # M_dot = dy/dx
        M_dot = np.array([ fit_list[i][0] for i in range(len(fit_list)) ])

        # tgrow = y_data[:-1]/M_dot[:-1]

        # x_data = x_data[:-1]/t_eddy[i]
        # x_data = x_data - x_data[0]

        tgrow = fit_y/M_dot

        fit_x = fit_x/t_eddy[i]
        fit_x = fit_x - fit_x[0]

        x_data = np.copy(fit_x)

        tgrow = np.abs(tgrow)

        # tgrow = sg.savgol_filter(tgrow, window_length=11, polyorder=3)

        # x_data = x_data[tgrow>1e-3]
        # tgrow  = tgrow [tgrow>1e-3]

        if MHD_flag:
            plt.plot(x_data,tgrow, color=line_col[i], linestyle=linestyle_list[i_MHD]\
            ,linewidth=4, label=r'R/l$_{\rm shatter} = $'+ f'{R_lsh[i]}')
        else:
            plt.plot(x_data,tgrow, color=line_col[i], linestyle=linestyle_list[i_MHD]\
            ,linewidth=4)


        if MHD_flag:
            plt.plot(x_data,tgrow, color='k', linestyle=linestyle_list[i_MHD]\
            ,linewidth=5,zorder=-2)

            plt.axhline(tg.tgrow(100, M, R_lsh[i], 40, ps.t_cool_cloud),\
                linestyle='dotted',color=line_col[i],linewidth=3,label=r'Theoretical t$_{\rm grow}$')

# plt.yscale('log')
plt.legend(loc='upper center')

# plt.xlim(0,3)
# plt.ylim(1e-1,4e0)
# plt.ylim(0,1)
plt.ylim(0,0.5)

plt.xlabel(r't/t$_{\rm eddy}$')
plt.ylabel(r't$_{\rm grow}$ = $\frac{m_{\rm cold}}{\dot{m}_{\rm cold}}$')