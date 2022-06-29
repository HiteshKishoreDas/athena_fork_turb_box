import numpy as np
from matplotlib import pyplot as plt

class hst_data:

    def __init__(self, hst_arr, ncells):

        self.time  = hst_arr[:,0]
        self.dt    = hst_arr[:,1]

        self.mass_tot  = hst_arr[:,2]

        self.mom1  = hst_arr[:,3]
        self.mom2  = hst_arr[:,4]
        self.mom3  = hst_arr[:,5]

        self.KE1   = hst_arr[:,6]
        self.KE2   = hst_arr[:,7]
        self.KE3   = hst_arr[:,8]
        self.E_tot = hst_arr[:,9]

        self.cold_gas = hst_arr[:,10]

        self.rho_avg    = hst_arr[:,11]/ncells
        self.rho_sq_avg = hst_arr[:,12]/ncells

        self.KE_tot = self.KE1+self.KE2+self.KE3

        self.turb_vel = np.sqrt(self.KE_tot*2/self.mass_tot)

        self.cs_avg    = hst_arr[:,13]/ncells
        self.tcool_avg = hst_arr[:,14]/ncells

        self.clumping_factor = self.rho_sq_avg/self.rho_avg**2

        self.cold_gas_fraction = self.cold_gas/self.mass_tot



if __name__=="__main__":
    
    hist = np.loadtxt("../work_turb_cloud/Turb.hst")
    # hist = np.loadtxt("../work_turb_cloud_Bfield/Turb.hst")
    # hist = np.loadtxt("Saved_videos/Turb.hst")
    # hist = np.loadtxt("../work_turb_dedt-mach/dedt_0.0001_box_1.0/Turb.hst")
    # hist = np.loadtxt("../work_turb_dedt-mach/dedt_0.0001_box_2.0/Turb.hst")


    cloud_time = 40


    time  = hist[:,0]
    dt    = hist[:,1]

    mass_tot  = hist[:,2]

    mom1  = hist[:,3]
    mom2  = hist[:,4]
    mom3  = hist[:,5]

    KE1   = hist[:,6]
    KE2   = hist[:,7]
    KE3   = hist[:,8]
    E_tot = hist[:,9]

    cold_gas = hist[:,10]

    rho_sum    = hist[:,11]
    rho_sq_sum = hist[:,12]


    KE_tot = KE1+KE2+KE3

    plt.figure(figsize=(8,8))
    plt.yscale("log")
    plt.plot(time,KE_tot)
    plt.xlabel("t")
    plt.ylabel("KE")

    plt.axvline(cloud_time,color='k',linestyle='dashed')

    turb_vel = np.sqrt(KE_tot*2/mass_tot)

    plt.figure(figsize=(8,8))
    # plt.yscale("log")
    plt.plot(time,turb_vel,\
        linewidth=3)
    plt.xlabel("t")
    plt.ylabel("Avg. v")
    plt.axvline(cloud_time,color='k',linestyle='dashed')

    N = 128
    dim = 3
    N_cells = N**dim

    clumping_factor = rho_sq_sum*N_cells/rho_sum

    plt.figure(figsize=(8,8))
    plt.yscale("log")
    plt.plot(time,clumping_factor,\
        linewidth=3)
    plt.xlabel("t")
    plt.ylabel("Clumping factor")
    plt.axvline(cloud_time,color='k',linestyle='dashed')

    plt.figure(figsize=(8,8))
    # plt.yscale("log")
    plt.plot(time, cold_gas,\
        linewidth=3)
    plt.xlabel("t")
    plt.ylabel("Cold gas")
    plt.axvline(cloud_time,color='k',linestyle='dashed')

    plt.figure(figsize=(8,8))
    # plt.yscale("log")
    plt.plot(time,cold_gas/mass_tot,\
        linewidth=3)
    plt.xlabel("t")
    plt.ylabel("Cold gas fraction")
    plt.axvline(cloud_time,color='k',linestyle='dashed')

    # N = 201
    # vel_smooth = np.convolve(turb_vel, np.ones(N)/N, mode='valid')

    # plt.plot(hist[int(N/2):-int(N/2),0],vel_smooth,linewidth=2)