import numpy as np
import v_turb as vt

# import sys
import os

import globals as g
import sys

cwd = os.getcwd()
repo_abs_path = cwd[: -len(cwd.split("/")[-1])]

cooling_dir = repo_abs_path + "cooling_scripts/"

sys.path.append(cooling_dir)
import cooling_fn as cf

parent_dir = cwd.split("/")[-1]

### FOR INPUT FILE

cloud_chi_temp = 100

cloud_pos_x = 0.0
cloud_pos_y = 0.0
cloud_pos_z = 0.0

# * Relevant temperature parameters

T_floor = 40000  # No cooling below T_floor
T_ceil = 100000000  # No gas above T_ceil

# T_hot     = T_hot_req*(0.7/2)*(4/7.5)*(4/7.5)*0.25   # Initial ambient hot medium temperature
T_hot_req = T_floor * cloud_chi_temp  # Required T_hot
T_hot = T_hot_req  # Initial ambient temperature

T_cold = 2 * T_floor  # for cold gas mass calculation

T_cut_mul = 0.5  # for cooling cutoff
T_cut = T_cut_mul * T_hot_req

amb_rho = 1.0  # ambient density

P_floor = 5 * 1e-4  # pressure floor

# * Chemical composition

Zsol = 1.0
Xsol = 1.0

X = Xsol * 0.7381
Z = Zsol * 0.0134
Y = 1 - X - Z

mu = 1.0 / (2.0 * X + 3.0 * (1.0 - X - Z) / 4.0 + Z / 2.0)
mue = 2.0 / (1.0 + X)
muH = 1.0 / X
mH = 1.0

# R_lsh = np.array([5,10,50,250,500])         # For M = 0.25
# R_lsh = np.array([10,50,100,250,500,1000])     # For M = 0.5
# R_lsh = np.array([250,500,1000,2500,5000])  # For M = 0.9

R_lsh = np.array([310])  # For M = 0.5

# * Cooling times
t_cool_cloud = cf.tcool_calc(amb_rho * cloud_chi_temp, T_floor, Zsol)
t_cool_mix = cf.tcool_calc(
    amb_rho * np.sqrt(cloud_chi_temp), np.sqrt(T_floor * T_hot_req), Zsol
)
t_cool_amb = cf.tcool_calc(amb_rho, T_hot_req, Zsol)
t_cool_cut = cf.tcool_calc(amb_rho / T_cut_mul, T_cut_mul * T_hot_req, Zsol)


l_sh = vt.cs_calc(T_floor, mu) * t_cool_cloud  # Cooling length scale

cloud_radius_temp = R_lsh * l_sh

L_box = cloud_radius_temp * 40  # Box size


# Cooling flag
cooling_flag = 1  # 1 for cooling and 0 for no cooling
# Cooling() is added(not added) to Source() depending on the flag


# Cloud flag
cloud_flag = 1  # 1 for a cloud and 0 for no cloud
# Cloud_init() is added(not added) to Source() depending on the flag

# Magnetic field flag
B_flag = 1  # 1 for adding magnetic fields

Mach_arr = np.array([0.25, 0.5, 0.9])

M = 0.5  # Required Mach number
i_mach = np.argwhere(Mach_arr == M)[0][0]


if B_flag:
    multi = 2
    dedt = vt.dedt_calc_hydro(M, amb_rho, T_hot_req, L_box) * multi

else:
    dedt = vt.dedt_calc_hydro(M, amb_rho, T_hot_req, L_box)

# Box sizes
x1max = L_box / 2
x1min = -1 * x1max

x2max = L_box / 2
x2min = -1 * x2max

x3max = L_box / 2
x3min = -1 * x3max

# Number of cells
nx1 = np.array([64])
nx2 = np.array([64])
nx3 = np.array([64])

nx1_mesh = np.array([16])
nx2_mesh = np.array([16])
nx3_mesh = np.array([16])

# predicted turbulent velocity
v_turb_predict = M * vt.cs_calc(T_hot_req, mu)

# t_eddy using the predicted turbulent velocity
t_eddy = L_box / v_turb_predict

tlim_trb = 7 * t_eddy
tlim_cld = 15 * t_eddy


# dt for history output from tlim
hst_dt_N = 2000

hst_dt_trb_arr = tlim_trb / hst_dt_N
hst_dt_cld_arr = tlim_cld / (hst_dt_N * (tlim_cld - tlim_trb) / tlim_trb).astype(int)

# dt for output files
out_dt_N = 500

out_dt_trb_arr = tlim_trb / out_dt_N
out_dt_cld_arr = tlim_cld / (out_dt_N * (tlim_cld - tlim_trb) / tlim_trb).astype(int)

# t_corr, t_corr ~ t_eddy
t_corr_arr = t_eddy

# dt_drive from t_eddy, using dt_drive << t_eddy
dt_drive_arr = t_eddy / 1000

# dt for restart files
rst_dt_N = 10

rst_dt_trb_arr = tlim_trb / rst_dt_N
rst_dt_cld_arr = tlim_cld / (rst_dt_N * (tlim_cld - tlim_trb) / tlim_trb).astype(int)

# ratio of shear component
f_shear = 1.0

# rseed
rseed = 1

# cloud properties
cloud_radius = cloud_radius_temp
cloud_chi = cloud_chi_temp
cloud_time = 7 * t_eddy


# Magnetic fields

# beta_list = np.array([1,2,5,10,100,1000])
# beta_list = np.array([2])

beta_list = 100

if B_flag:
    P_th = (T_hot_req / (g.KELVIN * g.mu)) * amb_rho
    P_B = P_th / beta_list

    B_mag = np.sqrt(P_B * 2.0)
    # Assuming that cgs relation between B_mag and P_B is used

    B_x = 0.0  # np.zeros_like(B_mag)
    B_y = 0.0  # np.zeros_like(B_mag)
    B_z = B_mag

else:
    B_x = 0.0
    B_y = 0.0
    B_z = 0.0

### FOR JOB SCRIPT

n_cores = (nx1 * nx2 * nx3) / (nx1_mesh * nx2_mesh * nx3_mesh)

# cluster_name = "cobra"
cluster_name = "freya"

if cluster_name == "freya":
    dir_path_add = "mpa/"
else:
    dir_path_add = ""

# queue = "p.24h" #"medium" # "n0064"
queue = "p.24h"
ntasks_per_node = 32

nodes = (n_cores / ntasks_per_node).astype(int)

time_limit_turb = ["04:49:00"]  # ["23:49:00"] #["03:00:00","12:00:00"]#,"23:49:00"]
time_limit_turb_rst = ["04:35:00"]  # ["23:35:00"] #["02:45:00","11:45:00"]#,"23:35:00"]

time_limit_cloud = ["04:59:00"]  # ["23:59:00"] #["03:00:00","12:00:00"]#,"23:49:00"]
time_limit_cloud_rst = [
    "04:55:00"
]  # ["23:55:00"] #["02:45:00","11:45:00"]#,"23:35:00"]

restart_N = 10

### TEST PARAMETERS

t_cc = np.sqrt(cloud_chi) * cloud_radius / v_turb_predict


### FILE NAME


def filename_turb_add(i, j):
    if B_flag:
        return f"_Rlsh{i}_{R_lsh[i]}_res{j}_{nx1[j]}_rseed_{rseed}_M_{M}_fshear_{f_shear}_beta_{beta_list}"
    else:
        return f"_Rlsh{i}_{R_lsh[i]}_res{j}_{nx1[j]}_rseed_{rseed}_M_{M}_fshear_{f_shear}_hydro"


def filename_cloud_add(i, j):
    if B_flag:
        return f"_Rlsh{i}_{R_lsh[i]}_res{j}_{nx1[j]}_rseed_{rseed}_M_{M}_chi_{cloud_chi}_beta_{beta_list}"
    else:
        return f"_Rlsh{i}_{R_lsh[i]}_res{j}_{nx1[j]}_rseed_{rseed}_M_{M}_chi_{cloud_chi}_hydro"


def filename_cloud_func(i, j, rseed, Mach, cloud_chi, beta, MHD_flag):
    if MHD_flag:
        return f"_Rlsh{i}_{R_lsh[i]}_res{j}_{nx1[j]}_rseed_{rseed}_M_{Mach}_chi_{cloud_chi}_beta_{beta}"
    else:
        return f"_Rlsh{i}_{R_lsh[i]}_res{j}_{nx1[j]}_rseed_{rseed}_M_{Mach}_chi_{cloud_chi}_hydro"
