import numpy as np
import sys
import para_scan as ps

sys.path.insert(0, '../vis/python')
import athena_read as ar
# sys.path.insert(0, '../../cooling_power_law')
# sys.path.insert(0, '../')



ncells = 64*64*640
# amb_rho = 1
# chi = 100
cold_rho = ps.amb_rho*ps.chi_cold

def lum_f (emissivity):
    return np.sum(emissivity)

def B_f (B, B0):
    return np.sqrt(np.sum((B - B0)**2)/ncells)

def entangle(B_vector, rho, i,j,k):

    #* Boolean matrix for mixed gas
    rho_mix = (rho>2*ps.amb_rho[i])*(rho<cold_rho[i]/2)

    #* Magnetic field in mixed gas
    B_mix_x = B_vector[0]*rho_mix
    B_mix_y = B_vector[1]*rho_mix
    B_mix_z = B_vector[2]*rho_mix

    #* Dot product of B_mix with neigbour in z direction
    B_dot = B_mix_x*np.roll(B_mix_x,-1,axis=2)+\
            B_mix_y*np.roll(B_mix_y,-1,axis=2)+\
            B_mix_z*np.roll(B_mix_z,-1,axis=2)

    return np.sum(np.sqrt(B_dot))
                 
def loop_func(i,j,k, B_flag):

    file_add = ps.filename_mix_add_ext(i,j,k,B_flag)
    dir_name = f'mix{file_add}'

    luminosity = []
    time_list = []
    B = []
    entanglement = []

    luminosity.append([0])
    time_list.append([0])
    B.append([0])
    entanglement.append([0])

    for n in range(0,200):

        print(f'n: {n}')

        file_name_lum =f'{dir_name}/Turb.out4.{str(n).zfill(5)}.athdf'

        try:
            data = ar.athdf(file_name_lum)
        except:
            break

        emissivity = -1*(data['user_out_var0'])

        luminosity.append(lum_f(emissivity))
        time_list.append(data['Time'])

        file_name_prim =f'{dir_name}/Turb.out2.{str(n).zfill(5)}.athdf'

        try:
            data = ar.athdf(file_name_prim)
        except:
            break

        if B_flag:
            B.append( [B_f(data['Bcc1'],data0['Bcc1']),\
                       B_f(data['Bcc2'],data0['Bcc2']),\
                       B_f(data['Bcc3'],data0['Bcc3'])]  )

            entanglement.append( entangle([data['Bcc1'], data['Bcc2'], data['Bcc3']], data['rho'] ,i,j,k ))

    return luminosity, B, time_list, entanglement


for B_fl in [True, False]:
    for i in range(len(ps.amb_rho)):
        for j in range(len(ps.Ma)):
            for k in range(3):

                if not(B_fl) and k!=0:
                    break

                file_add = ps.filename_mix_add_ext(i,j,k,B_fl)
                dir_name = f'mix{file_add}'

                print(f'Analysing {dir_name} ...')

                file0 = f'{dir_name}/Turb.out2.{str(0).zfill(5)}.athdf'
                data0 = ar.athdf(file0)

                luminosity, B, time_list, entanglement = loop_func(i,j,k, B_fl)

                luminosity = luminosity[1:]
                time_list  = time_list[1:]
                B = B[1:]
                entanglement = entanglement[1:]

                np.savetxt(f"save_arr/luminosity{file_add}", np.array(luminosity))
                np.savetxt(f"save_arr/time_list{file_add}", np.array(time_list))
                np.savetxt(f"save_arr/B{file_add}", np.array(B))
                np.savetxt(f"save_arr/entanglement{file_add}", np.array(entanglement))
