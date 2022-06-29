import numpy as np
import yt

from globals import *

def hdf2arr_fp(file_path, field="rho"):

    ds = yt.load(file_path)
    time = ds.current_time

    all_data_level_0 = ds.covering_grid(
        level=0, left_edge=[0, 0.0, 0.0], dims=ds.domain_dimensions
    )

    return np.array(all_data_level_0[field]), float(str(time).split()[0])
def hdf2arr_ds(ds, field="rho"):

    # ds = yt.load(file_path)
    time = ds.current_time

    all_data_level_0 = ds.covering_grid(
        level=0, left_edge=[0, 0.0, 0.0], dims=ds.domain_dimensions
    )

    return np.array(all_data_level_0[field]), float(str(time).split()[0])

if __name__=="__main__":

    from matplotlib import pyplot as plt

    CONST_pc  = 3.086e18
    CONST_yr  = 3.154e7
    CONST_amu = 1.66053886e-24
    CONST_kB  = 1.3806505e-16
    unit_length = CONST_pc*1e3  # 1 kpc
    unit_time   = CONST_yr*1e6  # 1 Myr
    unit_density = CONST_amu    # 1 mp/cm-3
    unit_velocity = unit_length/unit_time
    KELVIN = unit_velocity*unit_velocity*CONST_amu/CONST_kB

    g = 5/3
    T_floor = 10000.0
    T_ceil = 1e8

    X = 1.0
    Y = 0.0
    Z = 0.0

    mu = 1.0/(2.*X+ 3.*(1.-X-Z)/4.+ Z/2.)
    mue = 2.0/(1.0+X)
    muH = 1.0/X
    mH = 1.0

    # data_dir = "../work_turb_cloud/output/"
    data_dir = "../work_turb_cloud/"

    ds = yt.load(data_dir+"Turb.out2.00002.athdf")
    rho_arr, time = hdf2arr_fp(data_dir+"Turb.out2.00002.athdf")
    prs_arr, time = hdf2arr_fp(data_dir+"Turb.out2.00002.athdf",field="press")



    T_arr = (prs_arr/rho_arr)*KELVIN*mu

    gamma = 5/3

    cs_arr = np.sqrt(gamma*prs_arr/rho_arr)

    print("Average c_s in cloud  : ",np.average(cs_arr[rho_arr>20]))
    print("Average c_s in hot gas: ",np.average(cs_arr[rho_arr<5]))

    # print

    plt.figure(figsize=(10,10))
    plt.imshow(np.log10(rho_arr)[:,:,10])
    plt.colorbar()
    plt.show()

    print(f"Time: {time}")

    import h5py as hp

    df = hp.File(data_dir+"Turb.out2.00002.athdf",'r')

    prim = df['prim']

    print(prim.shape)