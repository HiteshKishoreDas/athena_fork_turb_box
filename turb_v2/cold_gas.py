import numpy as np
import matplotlib.pyplot as plt

import para_scan as ps
import athena_read as ar
from globals import *

for i in range(len(ps.R_lsh)):

    fn_suffix = ps.filename_cloud_func(i,j=0,rseed=1,Mach=0.5,cloud_chi=100,beta=100,MHD_flag=False)

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

        cold_gas.append(np.sum(data['rho'][T<1.5*ps.T_cold]))
        time.append(data['Time'])

    save_arr = []
    save_arr.append(time)
    save_arr.append(cold_gas)

    np.savetxt(f"save_arr/save_arr{fn_suffix}", np.array(save_arr))
