
#%%
import numpy as np
from matplotlib import markers, pyplot as plt
import yt

from globals import *

import hdf5_to_nparray as hn



def cold_mass(ds,T_cold = 1e5,dvol=1.0):
    T_arr, time = hn.hdf2arr_ds(ds,field="temp")
    rho_arr, time = hn.hdf2arr_ds(ds,field="density")
    
    cold_gas_t = np.sum(rho_arr[T_arr<T_cold])*dvol

    return cold_gas_t, float(str(time).split()[0])

def avg_temp(ds, mass_weighted_flag=False, dvol = 1.0):
    rho_arr, time = hn.hdf2arr_ds(ds,field="rho")
    T_arr, time = hn.hdf2arr_ds(ds,field="temp")

    if mass_weighted_flag:
        return np.average(T_arr,weights=rho_arr*dvol)
    else:
        return np.average(T_arr)

def total_mass(ds,dvol=1.0):
    rho_arr, time = hn.hdf2arr_ds(ds,field="rho")
    T_arr, time = hn.hdf2arr_ds(ds,field="temp")

    tot_gas_t = np.sum(rho_arr)*dvol

    return tot_gas_t, float(str(time).split()[0])

def cold_gas_fraction (ds,T_cold = 1e5,dvol=1.0):

    cold_gas,time = cold_mass(ds,T_cold,dvol)
    total_gas,time = total_mass(ds,dvol)

    return cold_gas/total_gas, float(str(time).split()[0])


#%%

if __name__=="__main__":

    data_dir = "../work_turb_cloud/"

    ts = yt.load(data_dir+"Turb.out2.*.athdf")

    cold_gas_frac = []
    time_list = []
    T_avg = []

    for ds in ts:

        cold_gas_frac_t,time = cold_gas_fraction(ds)

        cold_gas_frac.append(cold_gas_frac_t)
        time_list.append(time)

        T_avg.append(avg_temp(ds,mass_weighted_flag=True))


    plt.figure()
    # plt.yscale('log')

    plt.plot(range(len(cold_gas_frac)),cold_gas_frac)
    plt.scatter(range(len(cold_gas_frac)),cold_gas_frac)
    # plt.scatter(time_list,cold_gas_frac)
    # plt.xlim(2,5)
    plt.savefig("cold_gas/cold_gas_fraction.png")

    plt.figure()
    plt.yscale('log')
    plt.plot(range(len(T_avg)),T_avg)
    plt.scatter(range(len(T_avg)),T_avg)
    plt.savefig("cold_gas/average_temperature.png")

    # ts = yt.load(data_dir+"Turb.out2.*.athdf")

    # rho_arr, time = hn.hdf2arr_ds(ts[4],field="rho")
    # prs_arr, time = hn.hdf2arr_ds(ts[4],field="press")

    # T_arr = prs_arr/rho_arr

    # T_cold = 0.05

    # T_slice = T_arr[:,:,20]

    # T_slice[T_slice<=T_cold] = 0.0
    # T_slice[T_slice>T_cold] = 1.0

    # plt.figure()
    # plt.imshow(T_slice)
    # plt.colorbar()
    # plt.show()

    # print(f"Time: {time}")

    # plt.figure()

    # for i in range(np.shape(T_arr)[2]):
    #     T_slice_i = T_arr[:,:,i]
    #     plt.scatter(i,np.sum(T_slice_i<=T_cold),color='tab:blue')

    # plt.show()

# %%
