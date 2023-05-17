import numpy as np
import matplotlib.pyplot as plt

import para_scan as ps
import athena_read as ar
from globals import *

MHD_list = [True, False]

# M_list = [0.25, 0.5, 0.9]
M_list = [0.5]

for M in M_list:
    for mhd in MHD_list:
        # for i in range(len(ps.R_lsh)):
        for i in range(len(ps.R_lsh)): 

            # fn_suffix = ps.filename_cloud_func(i,j=0,rseed=1,Mach=M,cloud_chi=100,beta=100,MHD_flag=mhd)
            if mhd:
                fn_suffix = f'_Rlsh{i}_{ps.R_lsh[i]}_res0_256_rseed_1_M_{M}_chi_100_beta_100'
            else:
                fn_suffix = f'_Rlsh{i}_{ps.R_lsh[i]}_res0_256_rseed_1_M_{M}_chi_100_hydro'

            cold_gas = []
            time = []

            print(fn_suffix)

            for N in range(501,566):

                print(f'Snapshot: {N}')

                fn = f"para_scan{fn_suffix}/Turb.out2.00{N}.athdf"

                # print(fn)

                try:
                    data = ar.athdf(fn)
                except:
                    print("Last snapshot!")
                    break


                T = (data['press']/data['rho'])*KELVIN*mu

                cold_gas.append(np.sum(data['rho'][T<ps.T_cold]))
                time.append(data['Time'])

            save_arr = []
            save_arr.append(time)
            save_arr.append(cold_gas)

            np.savetxt(f"save_arr/save_arr{fn_suffix}", np.array(save_arr))
