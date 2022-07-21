#%%
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib as mt

plt.style.use('../plot_scripts/plot_style.mplstyle')

import v_turb as vt
import para_scan as ps


beta_list = np.array([100])
Rlsh_list = np.array([100,250,500,1000,2500, 50, 10, 5, 2, 1])
res_list  = np.array([128])

file_list  = ['mhd','hydro']
MHD_or_not = [True,False]
rseed_list = [1]
linestyle_list = ['solid','dashed']

marker_list_hyd = ['X','o','D','^','v']
marker_list_mhd = ['o','D','^','v','X']

gamma = 5/3

M_arr = [0.25,0.5,0.9]

survival = []
Mach_list = []
R_list = []
MHD_flag_list = []

for i_m, m in enumerate(M_arr):
    for i_fl,file_name in enumerate(file_list):
        for i_b, beta in enumerate(beta_list):
            for i_r, Rlsh in enumerate(Rlsh_list):
                for res in res_list:
                    for rseed in rseed_list:


                        v_turb_predict = m*vt.cs_calc(ps.T_hot_req,ps.mu)
                        t_eddy_m = Rlsh*ps.l_sh*40/v_turb_predict

                        fn_suffix = ps.filename_cloud_func(i_r,j=0,rseed=1,Mach=m,cloud_chi=100,beta=100,MHD_flag=MHD_or_not[i_fl])
                
                        try:
                            save_arr = np.loadtxt(f"save_arr/save_arr{fn_suffix}")
                            time     = save_arr[0,:]
                            cold_gas = save_arr[1,:]
                        except:
                            continue

                
                        x_data = time/t_eddy_m
                        y_data = cold_gas/cold_gas[0]

                        R_list.append(Rlsh)
                        Mach_list.append(m)
                        MHD_flag_list.append(MHD_or_not[i_fl])

                        # surv_cut = 2.0

                        survival.append(y_data[-1])

                        # if y_data[-1]>surv_cut:
                        #     survival.append(True)
                        # else:
                        #     survival.append(False)

survival = np.array(survival)

#%%

MHD_mark = ['D' if s==True else 'o' for s in MHD_flag_list]
MHD_size = [100  if s==True else 600 for s in MHD_flag_list]
MHD_zord = [1   if s==True else -1 for s in MHD_flag_list]

# surv_col = ['tab:blue' if s==True else 'tab:red' for s in survival]

#%%

colorbar_max =  0.3
colorbar_min = -0.3

# cmap_name = 'viridis'
cmap_name = 'bwr_r'
cmap = mt.cm.get_cmap(cmap_name)

cb_qnt = np.copy(np.log10(survival))
norm = plt.Normalize(vmin=colorbar_min,vmax=colorbar_max)

# line_col = cmap(np.log10(beta_list)/np.log10(beta_list).max())
line_col = cmap(norm(cb_qnt))

# Add labels to the plot



#%%


fig2, ax2 = plt.subplots(1,1,figsize=(10,10))

for i_s in range(len(survival)):

    if MHD_flag_list[i_s]:
        ax2.scatter(Mach_list[i_s], R_list[i_s], \
            edgecolor='w', label='MHD', linewidth=2,\
            color=line_col[i_s], marker='D', s=MHD_size[i_s], zorder=MHD_zord[i_s])

    else:
        ax2.scatter(Mach_list[i_s], R_list[i_s], \
            edgecolor='k', label='HD', linewidth=2,\
            color=line_col[i_s], marker='o', s=MHD_size[i_s], zorder=MHD_zord[i_s])

x_line = np.linspace(0,1.0,num=100)

mul1 = (4/3) * (ps.t_cool_mix/ps.t_cool_cloud)
mul2 = (ps.t_cool_mix/ps.t_cool_cloud)

ax2.plot(x_line, x_line*10**(x_line) * mul1,linestyle='dashed',color='k')
ax2.plot(x_line, x_line*3**(x_line) * mul2,linestyle='dotted',color='k')

style_2 = dict(size=20, color='k', rotation=0)
ax2.text(0.05, 100, \
    # r"$ \frac{R}{l_{\rm shatter}} = \frac{4}{3}~\frac{t_{\rm cool,mix}}{t_{\rm cool, floor}}~\mathcal{M}~10^{\mathcal{M}}$",\
    r"$ \frac{4}{3}~\frac{t_{\rm cool,mix}}{t_{\rm cool, floor}}~\mathcal{M}~10^{\mathcal{M}}$",\
    **style_2)
ax2.text(0.4, 20, \
    # r"$ \frac{R}{l_{\rm shatter}} = \frac{t_{\rm cool,mix}}{t_{\rm cool, floor}}~\mathcal{M}~3^{\mathcal{M}}$",\
    r"$ \frac{t_{\rm cool,mix}}{t_{\rm cool, floor}}~\mathcal{M}~3^{\mathcal{M}}$",\
    **style_2)

ax2.set_yscale('log')
ax2.set_ylim(1e-1,1e4)
ax2.set_xlim(0,)

ax2.set_xlabel(r'$\mathcal{M}$')
ax2.set_ylabel(r'$\frac{R_{\rm cl}}{l_{\rm shatter}}$', rotation=0)

ax2.set_title(f"Resolution: {res_list[0]}"+r'$^{3}$')

legend_elements = [Line2D([0], [0], marker='o', color='w', label='HD',
                          markerfacecolor=line_col[0], markersize=15, markeredgewidth=2, markeredgecolor='k'),
                #    Line2D([0], [0], marker='o', color='w', label='HD : Survived ',
                #           markerfacecolor='tab:blue', markersize=15, markeredgewidth=2, markeredgecolor='k'),
                   Line2D([0], [0], marker='D', color='w', label='MHD',
                          markerfacecolor=line_col[0], markersize=10, markeredgewidth=2, markeredgecolor='w')]
                #    Line2D([0], [0], marker='D', color='w', label='MHD: Survived ',
                #           markerfacecolor='tab:blue', markersize=10, markeredgewidth=2, markeredgecolor='w')]

plt.legend(handles=legend_elements, bbox_to_anchor=(0.99,0.2))
# plt.legend()
sm = plt.cm.ScalarMappable(cmap=cmap_name, norm=plt.Normalize(vmin=colorbar_min, vmax=colorbar_max))
plt.colorbar(sm, ax=ax2, label=r'log$_{10} \frac{M_{\rm cold}}{M_{\rm cold,0}}$')

#%%


fig1, ax1 = plt.subplots(1,1,figsize=(10,10))

t_cc = np.sqrt(ps.cloud_chi)*np.array(R_list)*ps.l_sh/(np.array(Mach_list)*vt.cs_calc(ps.T_hot_req,ps.mu))

tcool_mix_tcc = ps.t_cool_mix/t_cc

for i_s in range(len(survival)):

    if MHD_flag_list[i_s]:
        ax1.scatter(Mach_list[i_s], tcool_mix_tcc[i_s], \
            edgecolor='w', label='MHD', linewidth=2,\
            color=line_col[i_s], marker='D', s=MHD_size[i_s], zorder=MHD_zord[i_s])

    else:
        ax1.scatter(Mach_list[i_s], tcool_mix_tcc[i_s], \
            edgecolor='k', label='HD', linewidth=2,\
            color=line_col[i_s], marker='o', s=MHD_size[i_s], zorder=MHD_zord[i_s])

x_line = np.linspace(0,1.0,num=100)
# ax1.plot(x_line,(0.75)*(10**(-1*x_line)),linestyle='dashed',color='k')
ax1.plot(x_line,(1.0)*(10**(-np.log10(3)*x_line)),linestyle='dotted',color='k')
ax1.plot(x_line,(0.75)*(10**(-1*x_line)),linestyle='dashed',color='k')

style_1 = dict(size=20, color='k')
ax1.text(0.53, 0.14, \
    r"$ \frac{t_{\rm cool, mix}}{t_{\rm cc}} = \frac{3}{4}~10^{-\mathcal{M}}$",\
    **style_1, rotation=-20)
ax1.text(0.6, 0.48, \
    r"$ \frac{t_{\rm cool, mix}}{t_{\rm cc}} = 3^{-\mathcal{M}}$",\
    **style_1, rotation=-10)

ax1.set_yscale('log')
# plt.xscale('log')

ax1.set_xlabel(r'$\mathcal{M}$')
ax1.set_ylabel(r'$t_{\rm cool,mix}/t_{\rm cc}$')

ax1.set_xlim(0,1)
# ax1.set_ylim(1e-2,2e1)

ax1.set_title(f"Resolution: {res_list[0]}"+r'$^{3}$')

legend_elements = [Line2D([0], [0], marker='o', color='w', label='HD',
                          markerfacecolor=line_col[0], markersize=15, markeredgewidth=2, markeredgecolor='k'),
                #    Line2D([0], [0], marker='o', color='w', label='HD : Survived ',
                #           markerfacecolor='tab:blue', markersize=15, markeredgewidth=2, markeredgecolor='k'),
                   Line2D([0], [0], marker='D', color='w', label='MHD',
                          markerfacecolor=line_col[0], markersize=10, markeredgewidth=2, markeredgecolor='w')]
                #    Line2D([0], [0], marker='D', color='w', label='MHD: Survived ',
                #           markerfacecolor='tab:blue', markersize=10, markeredgewidth=2, markeredgecolor='w')]

plt.legend(handles=legend_elements, bbox_to_anchor=(0.99,0.2))

# sm = plt.cm.ScalarMappable(cmap=cmap_name, norm=plt.Normalize(vmin=cb_qnt.min(), vmax=cb_qnt.max()))
sm = plt.cm.ScalarMappable(cmap=cmap_name, norm=plt.Normalize(vmin=colorbar_min, vmax=colorbar_max))
plt.colorbar(sm, ax=ax1, label=r'log$_{10} \frac{M_{\rm cold}}{M_{\rm cold,0}}$')

# %%
