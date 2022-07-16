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

cmap_name = 'Paired'

cmap = mt.cm.get_cmap(cmap_name)

# cb_qnt = np.scopy(ps.dedt)
cb_qnt = np.copy(ps.R_lsh)

# line_col = cmap(np.log10(beta_list)/np.log10(beta_list).max())
line_col = cmap(cb_qnt/cb_qnt.max())


M = 0.5
v_turb_predict = M*vt.cs_calc(ps.T_hot_req,ps.mu)
t_eddy = ps.L_box/v_turb_predict

plt.figure()

for i_MHD, MHD_flag in enumerate(MHD):
    for i in range(3,len(ps.R_lsh)):

        fn_suffix = ps.filename_cloud_func(i,j=0,rseed=1,Mach=0.5,cloud_chi=100,beta=100,MHD_flag=MHD_flag)

        save_arr = np.loadtxt(f"save_arr/save_arr{fn_suffix}")

        time     = save_arr[0,:]
        cold_gas = save_arr[1,:]

        x_data = time
        y_data = cold_gas/cold_gas[0]

        
        dy = np.roll(y_data,-1) - y_data
        dx = np.roll(x_data,-1) - x_data

        M_dot = dy/dx

        tgrow = y_data[:-1]/M_dot[:-1]

        x_data = x_data[:-1]/t_eddy[i]
        x_data = x_data - x_data[0]

        tgrow = np.abs(tgrow)

        tgrow = sg.savgol_filter(tgrow, window_length=5, polyorder=3)

        x_data = x_data[tgrow>1e-3]
        tgrow  = tgrow [tgrow>1e-3]

        if MHD_flag:
            plt.plot(x_data,tgrow, color=line_col[i], linestyle=linestyle_list[i_MHD]\
            ,linewidth=4, label=r'R/l$_{\rm shatter} = $'+ f'{ps.R_lsh[i]}')
        else:
            plt.plot(x_data,tgrow, color=line_col[i], linestyle=linestyle_list[i_MHD]\
            ,linewidth=4)


        if MHD_flag:
            plt.plot(x_data,tgrow, color='k', linestyle=linestyle_list[i_MHD]\
            ,linewidth=5,zorder=-2)

            plt.axhline(tg.tgrow(100, M, ps.R_lsh[i], 40, ps.t_cool_cloud),\
                linestyle='dotted',color=line_col[i],linewidth=3,label=r'Theoretical t$_{\rm grow}$')

plt.yscale('log')
plt.legend(loc='lower right')

# plt.xlim(0,3)
# plt.ylim(1e-3,2e2)

plt.xlabel(r't/t$_{\rm eddy}$')
plt.ylabel(r't$_{\rm grow}$ = $\frac{m_{\rm cold}}{\dot{m}_{\rm cold}}$')