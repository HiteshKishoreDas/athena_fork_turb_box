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

    #* Copied cumulative cooling data
    tot_cool = np.copy(hst.total_cooling)

    dcool = np.roll(tot_cool, -1) - tot_cool
    dt    = np.roll(hst.time, -1) - hst.time

    dcool = dcool   [hst.cold_gas_fraction>0.1]
    dt    = dt      [hst.cold_gas_fraction>0.1]
    time  = hst.time[hst.cold_gas_fraction>0.1]

    box_full = np.argwhere(hst.cold_gas_fraction>0.998)

    if len(box_full)!=0:
        dcool = dcool   [:np.min(box_full)]
        dt    = dt      [:np.min(box_full)]
        time  = hst.time[:np.min(box_full)]

    box_volume = ps.box_length*ps.box_width*ps.box_width
    # luminosity = (dcool/dt)[1:-1]*box_volume/(ps.box_width**2)

    dz = (ps.box_length[0]/ps.nx3[0])

    luminosity = (dcool/dt)[1:-1] * dz
    time       = time[1:-1]

    return time, luminosity



# fig, ax = plt.subplots(nrows=2, ncols=2, figsize = (20,25))
fig, ax = plt.subplots(nrows=len(ps.amb_rho), ncols=2, figsize = (20,60))

Da_list       = []
L_avg_list    = []
Q_theory_list = []

for i_p in range(len(ps.amb_rho)):

    i = i_p #(len(ps.amb_rho)-1) - i_p

    for j in range(len(ps.Ma)):
        for k in range(3):
            for B_fl in [True, False]:
            # for B_fl in [False]:

                #* Non-magnetic cases
                if not(B_fl) and k!=0:
                    continue

                # For colors and linestyles
                if not(B_fl):
                    plot_i = 3
                else:
                    plot_i = k


                #* Read history file 
                file_add = ps.filename_mix_add_ext(i_p,j,k,B_fl)
                if not B_fl:
                    setup_name = file_add
                dir_name = f'mix{file_add}'

                hst = hr.hst_data(f'{dir_name}/Turb.hst', ncells, B_fl)     

                #* Damkohler number calculation
                u    = (ps.v_shear**0.8) * (ps.box_width[0]/ps.t_cool_cloud[i])**0.2

                tcool_tKH = np.round(ps.t_cool_mix[i]/ps.t_KH[0], 4)

                Da = (ps.box_width[0]/u)/(ps.t_cool_mix[i])


                time, luminosity = lum_fn(hst)

                ax[i,j].plot(time, luminosity, color=col_list[plot_i], label=legend_list[plot_i])

                if not B_fl:
                    
                    Da_list.append(Da)
                    
                    L_avg = np.average(luminosity[-500:])
                    
                    ax[i,j].axhline(L_avg, linestyle='dashed', color=col_list[plot_i],\
                        label=r'L$_{\rm avg}$'+ f' = {np.round(L_avg,3)}')


                    if Da>2 :     # Fast cooling
                        p_e = 1
                        u_e = 3/4
                        L_e = 1/4
                        t_e = -1/4
                    elif Da<=2 :   # Weak cooling
                        p_e = 1
                        u_e = 1/2
                        L_e = 1/2
                        t_e = -1/2 

                    Q0 = 1
                    # Q0 = 8.8e-8  # in erg cm^-2 s^-1
                    # Q0 = Q0 * unit_Q


                    P     = ps.amb_rho[i]*ps.T_hot/ps.mu    # in kB cm^-3 K
                    u    *= g.unit_velocity/1e5

                    Q_theory  = Q0 
                    Q_theory *= (P/160)**p_e
                    Q_theory *= (u/30)**u_e 
                    Q_theory *= (ps.box_width[0]*1000/100)**L_e
                    Q_theory *= (ps.t_cool_cloud[i]/0.03)**t_e

                    #* Q/Q_0
                    print(f'Q_theory = {Q_theory}, Da = {Da}')  
                    Q_theory_list.append(Q_theory) 

                    print(f'L_avg = {L_avg}')
                    L_avg_list.append(L_avg)


                    print('_______________________________________________')

        ax[i,j].set_yscale('log')
        ax[i,j].legend(loc='lower right')

        ax[i,j].set_xlabel('time (Myr)')
        ax[i,j].set_ylabel('Luminosity')

        # ax[i,j].set_ylim(0.99,1)

        ax[i,j].set_title(r"$\mathcal{M}_{\rm a}$"+f" = {ps.Ma[j]}; " + r"t$_{\rm cool}$/t$_{\rm KH}$" + f" = {tcool_tKH}" + f"; Da = {np.round(Da, 3)}")#+f"\n {setup_name}")


plt.show()

#%%

plt.figure()

plt.xscale('log')
plt.yscale('log')

Q_0 = 5.66e-8             # in code units

# unit_Q = g.unit_energy * (g.unit_length**-2) * (g.unit_time**-1)

sim_plot    = []
theory_plot = []

sim_plot    = np.array(L_avg_list) / Q_0 
# sim_plot    = sim_plot/sim_plot[-1]

theory_plot = np.array(Q_theory_list) #* 1e4
# theory_plot = theory_plot/theory_plot[-1]

plt.plot(Da_list, sim_plot, label='Simulation luminosity')
plt.plot(Da_list, theory_plot, label=r'Q/Q$_0$ from Theory $\times 10^4$')
# plt.plot(Da_list, theory_plot/sim_plot, label='Theory')

plt.legend()

plt.show()
# %%
