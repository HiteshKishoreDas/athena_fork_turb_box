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

# ncells = 64*64*640
# ncells = 256*256*2560
ncells = 128*128*1280

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


Da_list_line  = []
Da_list_point = []
L_avg_list    = []
Q_strong_list = []
Q_weak_list   = []

L_scaled_list = []
Q_strong_scaled_list = []
Q_weak_scaled_list   = []

for j in range(len(ps.Ma)):
    for i in range(len(ps.box_width)):


        #* Damkohler number calculation
        # v_turb in km/s
        u  =  50 * (ps.M**(4/5))
        u *=  ( ps.cs_cold*g.unit_velocity/(15*1e5) )**0.8 
        u *=  ( ps.t_cool_cloud[0]/0.03 )**-0.1

        # tcool_tKH = np.round(ps.t_cool_mix[i]/ps.t_KH[0], 4)

        t_turb = ps.box_width[i] / (u*1e5/g.unit_velocity) 
        # t_cool = ps.t_cool_mix[i]
        t_cool = ps.t_cool_Da[0]
        Da = t_turb/t_cool


        Da_cut = 2 

        #* For strong cooling line
        p_e_strong= 1
        u_e_strong= 3/4
        L_e_strong= 1/4
        t_e_strong= -1/4

        #* For weak cooling line
        p_e_weak  = 1
        u_e_weak  = 1/2
        L_e_weak  = 1/2
        t_e_weak  = -1/2 

        # Q0 = 1
        Q0 = 5.66e-8             # in code units

        P     = ps.amb_rho*ps.T_hot/ps.mu    # in kB cm^-3 K

        P_term_strong = (P/160)**p_e_strong
        P_term_weak   = (P/160)**p_e_weak

        u_term_strong = (u/30)**u_e_strong 
        u_term_weak   = (u/30)**u_e_weak

        L_term_strong = (ps.box_width[i]*1000/100)**L_e_strong
        L_term_weak   = (ps.box_width[i]*1000/100)**L_e_weak

        t_term_strong = (ps.t_cool_cloud[0]/0.03)**t_e_strong
        t_term_weak   = (ps.t_cool_cloud[0]/0.03)**t_e_weak

        Q_strong = Q0 
        Q_strong *= P_term_strong 
        Q_strong *= u_term_strong 
        Q_strong *= L_term_strong 
        Q_strong *= t_term_strong

        Q_weak   = Q0 
        Q_weak   *= P_term_weak
        Q_weak   *= u_term_weak
        Q_weak   *= L_term_weak
        Q_weak   *= t_term_weak

        #* Q/Q_0
        Q_strong_list.append(Q_strong) 
        Q_weak_list.append(Q_weak) 

        Q_strong_scaled_list.append(Q_strong)#/P_term_strong/u_term_strong/L_term_strong/t_term_strong)
        Q_weak_scaled_list.append(Q_weak)#/P_term_weak/u_term_weak/L_term_weak/t_term_weak)

        print(f'Da: {Da}')
        Da_list_line.append(Da)

        #* Read history file 
        file_add = ps.filename_mix_add_ext(i,j,0,False)#[:-7]
        dir_name = f'mix{file_add}'
        
        # hst = hr.hst_data(f'test_1/{dir_name}/Turb.hst', ncells, False)
        try:
            hst = hr.hst_data(f'{dir_name}/Turb.hst', ncells, False)
            print(f'{dir_name} read ...')
        except:
            print(f'Could not read the hst file in {dir_name}! ... ')
            continue

        Da_list_point.append(Da)

        time, luminosity = lum_fn(hst, i)

        # L_avg = np.average(luminosity[-250:])
        L_avg = np.average(luminosity[-50:])

        print(f'L_avg: {L_avg}')
        L_avg_list.append(L_avg)

        constant = 3

        if Da<2:
            L_scaled_list.append(constant * L_avg)#/P_term_weak/u_term_weak/L_term_weak/t_term_weak)
        else:
            L_scaled_list.append(constant * L_avg)#/P_term_strong/u_term_strong/L_term_strong/t_term_strong/1e4)


plt.figure()

plt.xscale('log')
plt.yscale('log')


# unit_Q = g.unit_energy * (g.unit_length**-2) * (g.unit_time**-1)

sim_plot    = []
theory_plot = []

sim_plot    = np.array(L_avg_list)
# sim_plot    = sim_plot/sim_plot[-1]

theory_strong_plot = np.array(Q_strong_list)
theory_weak_plot   = np.array(Q_weak_list)

# theory_plot = theory_plot/theory_plot[-1]

# plt.plot(Da_list, sim_plot, label='Simulation luminosity')
# plt.plot(Da_list, theory_plot, label=r'Q/Q$_0$ from Theory $\times 10^4$')
# plt.plot(Da_list, theory_plot/sim_plot, label='Theory')

plt.scatter(Da_list_point, L_scaled_list, label='Simulation luminosity')
plt.plot(Da_list_line, Q_strong_scaled_list, linestyle='dashed', label=r'Q: strong cooling')
plt.plot(Da_list_line, Q_weak_scaled_list,   linestyle='dashed', label=r'Q: weak cooling')

# Da_arr = np.array(Da_list)
# Q0_fit = (5/3)*1e-8*  (Da_arr**-(1/3))

# plt.plot(Da_arr, Q0_fit)

plt.legend()

plt.show()
# plt.savefig('hydro_check.png')