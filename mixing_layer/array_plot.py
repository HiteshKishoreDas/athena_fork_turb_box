import numpy as np
import matplotlib.pyplot as plt

plt.style.use('../plot_scripts/plot_style.mplstyle')

import sys
sys.path.insert(0, '../vis/python')
import athena_read as ar
import globals as g
import para_scan as ps


N = 10 

i = 0
j = 0
k = 1
B_fl = True

file_add = ps.filename_mix_add_ext(i, j, k, B_fl)
dir_name = f'mix{file_add}'
# file_name = f'{dir_name}/Turb.out2.{str(N).zfill(5)}.athdf'
file_name = f'{dir_name}/Turb.out4.{str(N).zfill(5)}.athdf'

ds = ar.athdf(file_name)

# rho_mix = np.sqrt(ps.chi_cold)*ps.amb_rho[i]

# rho_plot = ds['rho']
# rho_plot[rho_plot<rho_mix/4] =0 
# rho_plot[rho_plot>rho_mix*4] =0

# rho_plot = np.log10(rho_plot)

# zoom_ia = 100
# zoom_ib = 300

# plt.figure()
# # plt.imshow(rho[:,32,:], vmin = ps.amb_rho[i], vmax=100*ps.amb_rho[i], cmap='plasma_r')
# plt.imshow(rho_plot[zoom_ia:zoom_ib,32,:], vmin = np.log10(rho_mix/2), vmax= np.log10(rho_mix*2), cmap='plasma_r')
# plt.title(f"rho: {ps.amb_rho[i]}, Ma: {ps.Ma[j]}")
# plt.colorbar()


# v_dot_B  = ds['vel1']*ds['Bcc1']
# v_dot_B += ds['vel2']*ds['Bcc2']
# v_dot_B += ds['vel3']*ds['Bcc3']

# v_mag = ds['vel1']**2 + ds['vel2']**2 + ds['vel3']**2
# B_mag = ds['Bcc1']**2 + ds['Bcc2']**2 + ds['Bcc3']**2

# v_mag = np.sqrt(v_mag)
# B_mag = np.sqrt(B_mag)

# vB_cos = v_dot_B/(v_mag*B_mag)

# rho = ds['rho']
# va = B_mag/np.sqrt(rho)

# Ma = v_mag/va

# # plt.figure()
# # plt.imshow(Ma[:,32,:], vmin=0, vmax=4, cmap="RdGy")
# # plt.imshow(Ma[100:300,32,:], vmin=0, vmax=4, cmap="RdGy")
# plt.colorbar()

# plt.figure()
# plt.imshow(v_dot_B[zoom_ia:zoom_ib,32,:], vmin=-0.005, vmax=0.005, cmap="RdGy")
# plt.colorbar()


# plt.figure()
# plt.imshow(vB_cos[zoom_ia:zoom_ib,32,:], vmin=-1, vmax=1, cmap="RdGy")
# plt.colorbar()

# plt.figure(figsize=(10,20))
# plt.pcolormesh(ds['x1f'],ds['x3f'],ds['rho'][:,32,:], vmin = ps.amb_rho[i], vmax=100*ps.amb_rho[i])
# plt.axis('scaled')
# plt.colorbar()

plt.figure(figsize=(10,20))
plt.pcolormesh(ds['x1f'],ds['x3f'],np.log10(ds['user_out_var0'][:,32,:]))
plt.axis('scaled')
plt.colorbar()