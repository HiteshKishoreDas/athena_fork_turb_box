import numpy as np
import matplotlib.pyplot as plt

plt.style.use('../plot_scripts/plot_style.mplstyle')

import cmasher as cmr


import sys
import os
sys.path.insert(0, '../vis/python')
import athena_read as ar
import globals as g
import para_scan as ps

def B_analysis(ds):
    v_dot_B  = ds['vel1']*ds['Bcc1']
    v_dot_B += ds['vel2']*ds['Bcc2']
    v_dot_B += ds['vel3']*ds['Bcc3']

    v_mag = ds['vel1']**2 + ds['vel2']**2 + ds['vel3']**2
    B_mag = ds['Bcc1']**2 + ds['Bcc2']**2 + ds['Bcc3']**2

    v_mag = np.sqrt(v_mag)
    B_mag = np.sqrt(B_mag)

    vB_cos = v_dot_B/(v_mag*B_mag)

    rho = ds['rho']
    va  = B_mag/np.sqrt(rho)
    
    Ma  = v_mag/va

    return vB_cos, v_dot_B, Ma 


cmap = cmr.bubblegum

# N = 10 

i = 0
j = 0 
k = 1
B_fl = True

for i in range(len(ps.Lambda_fac)):
    for j in range(len(ps.Ma)):
        for k in range(3):
            # for B_fl in [True, False]:
            for B_fl in [False]:
                
                if not(B_fl) and k!=0:
                   continue

                file_add = ps.filename_mix_add_ext(i, j, k, B_fl)
                dir_name = f'mix{file_add}'

                dir_create_list = ['','/luminosity','/T','/beta','/Bx','/By','/Bz']

                for d_create in dir_create_list:
                    os.system(f"mkdir Plots/{dir_name}{d_create}")



                for N in range(101):

                    print(N)


                    file_name_prm = f'{dir_name}/Turb.out2.{str(N).zfill(5)}.athdf'
                    file_name_lum = f'{dir_name}/Turb.out4.{str(N).zfill(5)}.athdf'


                    try:
                        ds_prm = ar.athdf(file_name_prm)
                        ds_lum = ar.athdf(file_name_lum)
                    except:
                        print("Breaking to next file...")
                        break


                    #* Luminosity
                    plt.figure(figsize=(10,20))
                    plt.pcolormesh(ds_lum['x1f'],ds_lum['x3f'],np.average(ds_lum['user_out_var0'], axis=1) , cmap=cmap)
                    plt.title(dir_name)
                    plt.axis('scaled')
                    plt.colorbar()
                    plt.savefig(f"Plots/{dir_name}/luminosity/luminosity_{str(N).zfill(5)}")
                    plt.close()

                    T = (ds_prm['press']/ds_prm['rho']) * g.KELVIN * g.mu
                    T = np.log10(T)

                    print(f"Min pressure: {np.min(ds_prm['press'])}")

                    #* Temperature
                    plt.figure(figsize=(10,20))
                    # plt.pcolormesh(ds_prm['x1f'],ds_prm['x3f'],T[:,32,:], vmin=np.log10(4e4), vmax = np.log10(4e6))
                    plt.pcolormesh(ds_prm['x1f'],ds_prm['x3f'],np.average(T, axis=1), vmin=np.log10(4e4), vmax = np.log10(4e6))
                    plt.title(dir_name)
                    plt.axis('scaled')
                    plt.colorbar()
                    plt.savefig(f"Plots/{dir_name}/T/T_{str(N).zfill(5)}")
                    plt.close()



                    if B_fl:

                        B_mag = ds_prm['Bcc1']**2 + ds_prm['Bcc2']**2 + ds_prm['Bcc3']**2
                        B_mag = np.sqrt(B_mag)

                        beta = ds_prm['press']/ (0.5* B_mag**2)

                        #* Beta plot
                        plt.figure(figsize=(10,20))
                        plt.pcolormesh(ds_prm['x1f'],ds_prm['x3f'],beta[:,32,:])
                        plt.title(dir_name)
                        plt.axis('scaled')
                        plt.colorbar()
                        plt.savefig(f"Plots/{dir_name}/beta/beta_{str(N).zfill(5)}.png")
                        plt.close()

                        #* Bcc1
                        plt.figure(figsize=(10,20))
                        plt.pcolormesh(ds_prm['x1f'],ds_prm['x3f'],ds_prm['Bcc1'][:,32,:])
                        plt.title(dir_name)
                        plt.axis('scaled')
                        plt.colorbar()
                        plt.savefig(f"Plots/{dir_name}/Bx/Bx_{str(N).zfill(5)}.png")
                        plt.close()

                        #* Bcc2
                        plt.figure(figsize=(10,20))
                        plt.pcolormesh(ds_prm['x1f'],ds_prm['x3f'],ds_prm['Bcc2'][:,32,:])
                        plt.title(dir_name)
                        plt.axis('scaled')
                        plt.colorbar()
                        plt.savefig(f"Plots/{dir_name}/By/By_{str(N).zfill(5)}.png")
                        plt.close()

                        #* Bcc3
                        plt.figure(figsize=(10,20))
                        plt.pcolormesh(ds_prm['x1f'],ds_prm['x3f'],ds_prm['Bcc3'][:,32,:])
                        plt.title(dir_name)
                        plt.axis('scaled')
                        plt.colorbar()
                        plt.savefig(f"Plots/{dir_name}/Bz/Bz_{str(N).zfill(5)}.png")
                        plt.close()


                for d_c in dir_create_list[1:]:

                    if B_fl:
                        video_command = f'ffmpeg -framerate 5 -i Plots/{dir_name}{d_c}{d_c}_%05d.png -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2, fps=25, format=yuv420p" Plots/{dir_name}{d_c}{d_c}_m{i}{j}{k}.mp4'
                    else:
                        if d_c not in ['/beta','/Bx','/By','/Bz']:
                            video_command = f'ffmpeg -framerate 5 -i Plots/{dir_name}{d_c}{d_c}_%05d.png -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2, fps=25, format=yuv420p" Plots/{dir_name}{d_c}{d_c}_h{i}{j}{k}.mp4'
                        else:
                            continue

                    os.system(video_command)
 