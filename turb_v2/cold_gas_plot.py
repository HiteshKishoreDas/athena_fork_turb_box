import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mt 
import scipy.signal as sg

plt.style.use('../plot_scripts/plot_style.mplstyle')

import para_scan as ps
import v_turb as vt
import athena_read as ar
from globals import *


def add_legend(plot_args_lst, legend_loc = 'lower left', default_plot_args = {'color' : 'black'},
               **kwargs):
    """Adds another legend to plot.

    Keyword Arguments:
    plot_args_lst      -- List of plotting arguments to show.
    legend_loc         -- Location of new legend (default 'best')
    default_plot_args  -- Arguments will be used for every item in new legend.
    kwargs             -- Will be passed to `plt.legend` of new legend.

    Example:
           > maxpy.plot.add_legend([{'ls' : '-', 'label' : '8'}, {'ls' : '--', 'label' : '16'}],
           >                   default_plot_args = {'c' : 'k'},
           >                   title = r'$l_{\rm cell} / r_{\rm cl}$')
    Will add a legend with two different lines (both black).
    """
    ax = plt.gca()
    leg = ax.get_legend()

    linelst = []
    for cargs in plot_args_lst:
        for k, v in default_plot_args.items():
            if k not in cargs:
                cargs[k] = v
        l, = plt.plot(np.nan, **cargs)
        linelst.append(l)

    o = kwargs.copy()
    if 'loc' not in o:
        o['loc'] = legend_loc
    if legend_loc == 'above':
        o['loc'] = 'lower left'
        o['bbox_to_anchor'] = (0.5, 1.01)

    plt.legend(handles = linelst, **o)
    if leg is not None:
        ax.add_artist(leg) # Add old legend


# R_lsh = np.array([5,10,50,250,500])         # For M = 0.25
R_lsh = np.array([10,50,100,250,500,1000])     # For M = 0.5
# R_lsh = np.array([250,500,1000,2500,5000])  # For M = 0.9

MHD = [True, False]
linestyle_list = ['solid', 'dashed']

cmap_name = 'Paired'

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

fig1, ax1 = plt.subplots(1,1,figsize=(10,10))

for i_MHD, MHD_flag in enumerate(MHD):
    for i in range(len(R_lsh)):
    # for i in [9,8,7,6,5,1,2,3]:#range(len(ps.R_lsh)):
    # for i in [4]:#range(len(ps.R_lsh)):
    # for i in range(1,len(ps.R_lsh)-1):

        # fn_suffix = ps.filename_cloud_func(i,j=0,rseed=1,Mach=M,cloud_chi=100,beta=100,MHD_flag=MHD_flag)

        if MHD_flag:
            fn_suffix = f'_Rlsh{i}_{R_lsh[i]}_res0_256_rseed_1_M_{M}_chi_100_beta_100'
        else:
            fn_suffix = f'_Rlsh{i}_{R_lsh[i]}_res0_256_rseed_1_M_{M}_chi_100_hydro'

        try:
            save_arr = np.loadtxt(f"save_arr/save_arr{fn_suffix}")
            time     = save_arr[0,:]
            cold_gas = save_arr[1,:]
        except:
            continue

        x_data = time/t_eddy[i]
        y_data = cold_gas/cold_gas[0]

        y_cut = 1e-4

        if cold_gas.min() <= y_cut:
            x_data -= x_data[0]
            y_data = sg.savgol_filter(y_data, window_length=11, polyorder=3)

            ind_cut = np.where(y_data<y_cut)[0][0]

            x_data = x_data[:ind_cut+1]
            y_data = y_data[:ind_cut+1]
            y_data[ind_cut] = 1e-5
        
        else:
            x_data -= x_data[0]
            y_data = sg.savgol_filter(y_data, window_length=11, polyorder=3)

        if MHD_flag:
            ax1.plot(x_data,y_data, color=line_col[i], linestyle=linestyle_list[i_MHD]\
            ,linewidth=4, label=r'R/l$_{\rm shatter} = $'+ f'{R_lsh[i]}')
        else:
            ax1.plot(x_data,y_data, color=line_col[i], linestyle=linestyle_list[i_MHD]\
            ,linewidth=4)


        if MHD_flag:
            ax1.plot(x_data,y_data, color='k', linestyle=linestyle_list[i_MHD]\
            ,linewidth=5,zorder=-2)


ax1.axhline(1.0, linestyle='dotted', color='k')

ax1.set_yscale('log')
ax1.legend(loc='upper left')

# plt.xlim(0,1)
# plt.ylim(1e-1,4e2)
ax1.set_ylim(4e-1,10)
# ax1.set_ylim(9e-1,1.1)
# ax1.set_xlim(0,0.25)

ax1.set_xlabel(r't/t$_{\rm eddy}$')
ax1.set_ylabel(r'M$_{\rm cold}$/M$_{\rm cold,initial}$')

ax1.set_title(r"$\mathcal{M} = $"+f"{M}, Resolution: {res}"+ r"$^{3}$")

add_legend([{'ls' : '-', 'label' : 'MHD'}, {'ls' : '--', 'label' : 'HD'}],
                 default_plot_args = {'c' : 'k'})

sm = plt.cm.ScalarMappable(cmap=cmap_name, norm=plt.Normalize(vmin=cb_qnt.min(), vmax=cb_qnt.max()))
plt.colorbar(sm, ax=ax1, label=r'log $\frac{R}{l_{\rm shatter}}$')