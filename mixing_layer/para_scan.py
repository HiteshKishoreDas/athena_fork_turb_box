import numpy as np
import v_turb as vt

# import sys
import os

import globals as g
import sys

cwd = os.getcwd()
repo_abs_path = cwd[:-len(cwd.split('/')[-1])]

cooling_dir = repo_abs_path+'cooling_power_law/'

sys.path.append(cooling_dir)
import cooling_fn as cf

parent_dir = cwd.split('/')[-1]


chi_cold = 100

# Relevant temperature
T_floor   = 40000                       # No cooling below T_floor
T_ceil    = 100000000                   # No gas above T_ceil
T_hot     = T_floor*chi_cold            # Required T_hot
T_cold    = 2*T_floor                   # For cold gas mass calculation
T_cut_mul = 0.5                         # For cooling cutoff
T_cut     = T_cut_mul*T_hot             # For cooling cutoff

amb_rho   = np.array([1.0,0.01])                           # Density of the ambient medium

# Pressure floor
P_floor = 5*1e-4

# Chemical composition
Xsol = 1.0;
Zsol = 1.0;

X = Xsol * 0.7381;
Z = Zsol * 0.0134;
Y = 1 - X - Z;

mu  = 1.0/(2.*X+ 3.*(1.-X-Z)/4.+ Z/2.)
mue = 2.0/(1.0+X)
muH = 1.0/X
mH  = 1.0


# R_lsh = np.array([100,250,500,1000,2500])
R_lsh = np.array([2500])

t_cool_cloud = cf.tcool_calc(amb_rho*chi_cold,T_floor,Z)
t_cool_mix   = cf.tcool_calc(amb_rho*np.sqrt(chi_cold),np.sqrt(T_floor*T_hot),Z)
t_cool_amb   = cf.tcool_calc(amb_rho,T_hot,Z)
t_cool_cut   = cf.tcool_calc(amb_rho/T_cut_mul,T_cut_mul*T_hot,Z)


amb_rho_fix = 1.0
t_cool_cloud_fix = cf.tcool_calc(amb_rho_fix*chi_cold,T_floor,Z)
l_sh = vt.cs_calc(T_floor,mu)*t_cool_cloud_fix

cloud_radius = R_lsh*l_sh

box_width  = cloud_radius*2
box_length = box_width*10

# Cooling flag
cooling_flag = 1  # 1 for cooling and 0 for no cooling
                  # Cooling() is added(not added) to Source() depending on the flag 


# Cloud flag
cloud_flag   = 0  # 1 for a cloud and 0 for no cloud
                  # Cloud_init() is added(not added) to Source() depending on the flag 

# Magnetic field flag
B_flag       = 1  # 1 for adding magnetic fields

M  = 0.5     # Required Mach number
Ma = np.array([0.1, 10])

# Box sizes
x1max = box_width
x1min = np.array([0.0])

x2max = box_width
x2min = np.array([0.0])

x3max = box_length*0.5  # 0.8
x3min = x3max - box_length

# Number of cells
nx1 = np.array([64 ])
nx2 = np.array([64 ])
nx3 = np.array([640])

nx1_mesh = np.array([32])
nx2_mesh = np.array([32])
nx3_mesh = np.array([32])

# Initial profile settings 

front_thickness = box_width/40
v_shear         = M*vt.cs_calc(T_hot,mu)

knx_KH = 1.0
kny_KH = 1.0
amp_KH = 0.01

t_KH = np.sqrt(chi_cold)*box_width/v_shear
tlim  =  10*t_KH  #np.max(np.array([5*t_eddy, 2*t_cool_amb]) )

# dt for history output from tlim
hst_dt_arr = tlim/2000

# dt for output files
out_dt_arr = tlim/200

# dt for restart files
rst_dt_arr = tlim/10

# rseed
rseed = 1 


# Magnetic fields

beta_list = (2/g.g) * (Ma/M)**2

if B_flag:
    B_dir     = ['x','y','z']
else:
    B_dir     = ['h']

B_x = np.array([ [ [0.0] *len(B_dir) ] *len(Ma)] *len(amb_rho) )
B_y = np.array([ [ [0.0] *len(B_dir) ] *len(Ma)] *len(amb_rho) )
B_z = np.array([ [ [0.0] *len(B_dir) ] *len(Ma)] *len(amb_rho) )

if B_flag:

    for i_r, rho in enumerate(amb_rho): 

        for i_b, beta in enumerate(beta_list):

            P_th  = (T_hot/(g.KELVIN*g.mu))*rho
            P_B = P_th/beta

            # Assuming that cgs relation between B_mag and P_B is used
            B_mag = np.sqrt(P_B * 2.0)

            for i_d, B_d in enumerate(B_dir): 
               if B_d=='x':
                   B_x[i_r, i_b, i_d] = B_mag
               elif B_d=='y':
                   B_y[i_r, i_b, i_d] = B_mag
               elif B_d=='z':
                   B_z[i_r, i_b, i_d] = B_mag
               else:
                   print('! Invalid magnetic field direction (B_dir)! ... ')
                   exit


### FOR JOB SCRIPT

n_cores = (nx1*nx2*nx3)/(nx1_mesh*nx2_mesh*nx3_mesh)

queue = "p.24h"
ntasks_per_node = 32

nodes = (n_cores/ntasks_per_node).astype(int)

time_limit     = ["23:55:00"] #["03:00:00","12:00:00"]#,"23:49:00"]
time_limit_rst = ["23:52:00"] #["02:45:00","11:45:00"]#,"23:35:00"]



### FILE NAME

def filename_mix_add (i,j,k):

    if B_flag:
        return f'_rho{i}_Ma{j}_B{k}'
    else:
        return f'_rho{i}_Ma{j}_Bnot_hydro'
    # return f'_res_256_Rlsh_0'

#* For access to plotting scripts
def filename_mix_add_ext (i, j, k, MHD_flag):

    if MHD_flag:
        return f'_rho{i}_Ma{j}_B{k}'
    else:
        return f'_rho{i}_Ma{j}_Bnot_hydro'
