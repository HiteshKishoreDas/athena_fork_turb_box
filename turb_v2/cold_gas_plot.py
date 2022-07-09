import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mt 
import scipy.signal as sg

plt.style.use('../plot_scripts/plot_style.mplstyle')

import para_scan as ps
import v_turb as vt
import athena_read as ar
from globals import *


def add_legend(plot_args_lst, legend_loc = 'best', default_plot_args = {'color' : 'black'},
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
        o['loc'] = 'lower center'
        o['bbox_to_anchor'] = (0.5, 1.01)

    plt.legend(handles = linelst, **o)
    if leg is not None:
        ax.add_artist(leg) # Add old legend



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
    for i in range(len(ps.R_lsh)):

        fn_suffix = ps.filename_cloud_func(i,j=0,rseed=1,Mach=0.5,cloud_chi=100,beta=100,MHD_flag=MHD_flag)

        save_arr = np.loadtxt(f"save_arr/save_arr{fn_suffix}")

        time     = save_arr[0,:]
        cold_gas = save_arr[1,:]

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
            plt.plot(x_data,y_data, color=line_col[i], linestyle=linestyle_list[i_MHD]\
            ,linewidth=4, label=r'R/l$_{\rm shatter} = $'+ f'{ps.R_lsh[i]}')
        else:
            plt.plot(x_data,y_data, color=line_col[i], linestyle=linestyle_list[i_MHD]\
            ,linewidth=4)


        if MHD_flag:
            plt.plot(x_data,y_data, color='k', linestyle=linestyle_list[i_MHD]\
            ,linewidth=5,zorder=-2)


plt.yscale('log')
plt.legend(loc='lower right')

plt.xlim(0,3)
plt.ylim(1e-3,2e2)

plt.xlabel(r't/t$_{\rm eddy}$')
plt.ylabel(r'M$_{\rm cold}$/M$_{\rm cold,initial}$')



add_legend([{'ls' : '-', 'label' : 'MHD'}, {'ls' : '--', 'label' : 'HD'}],
                 default_plot_args = {'c' : 'k'})