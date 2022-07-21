import sys
sys.path.insert(0, '../vis/python')
import athena_read
import numpy as np

# dir_name = ""
dir_name = ["", "low_beta/", "low_beta/", "low_beta/"]

file_list = ['test_neun_hyd_cool','test_sechs_Bx_cool','test_seben_By_cool','test_acht_Bz_cool']
# file_list = ['test_sechs_Bx_cool']#,'test_seben_By_cool','test_acht_Bz_cool']


ncells = 64*64*640
amb_rho = 1
chi = 100
cold_rho = amb_rho*chi

def lum_f (emissivity):
    return np.sum(emissivity)

def B_f (B, B0):
    return np.sqrt(np.sum((B - B0)**2)/ncells)

def entangle(B_vector, rho):

    #* Boolean matrix for mixed gas
    rho_mix = (rho>2*amb_rho)*(rho<cold_rho/2)

    #* Magnetic field in mixed gas
    B_mix_x = B_vector[0]*rho_mix
    B_mix_y = B_vector[1]*rho_mix
    B_mix_z = B_vector[2]*rho_mix

    #* Dot product of B_mix with neigbour in z direction
    B_dot = B_mix_x*np.roll(B_mix_x,-1,axis=2)+\
            B_mix_y*np.roll(B_mix_y,-1,axis=2)+\
            B_mix_z*np.roll(B_mix_z,-1,axis=2)

    return np.sum(np.sqrt(B_dot))
                 
def loop_func(fl, i_fl):


    luminosity = []
    time_list = []
    B = []
    entanglement = []

    luminosity.append([0])
    time_list.append([0])
    B.append([0])
    entanglement.append([0])

    for n in range(0,60):

        print(f'n: {n}')

        file_name = fl + f'/{dir_name[i_fl]}Turb.out4.{str(n).zfill(5)}.athdf'

        try:
            data = athena_read.athdf(file_name)
        except:
            break

        # emissivity = np.abs(data['user_out_var0'])

        emissivity = -1*(data['user_out_var0'])

        luminosity.append(lum_f(emissivity))
        time_list.append(data['Time'])

        file_name = fl + f'/{dir_name[i_fl]}Turb.out2.{str(n).zfill(5)}.athdf'

        try:
            data = athena_read.athdf(file_name)
        except:
            break

        if i_fl!=0:
            B.append( [B_f(data['Bcc1'],data0['Bcc1']),\
                       B_f(data['Bcc2'],data0['Bcc2']),\
                       B_f(data['Bcc3'],data0['Bcc3'])]  )

            entanglement.append( entangle([data['Bcc1'], data['Bcc2'], data['Bcc3']], data['rho']  ))

    return luminosity, B, time_list, entanglement


for i_fl,fl in enumerate(file_list):

    save_arr_dir = "low_beta/"

    print(f'i_fl: {i_fl}')

    file0 = fl + f'/{dir_name[i_fl]}Turb.out2.{str(0).zfill(5)}.athdf'
    data0 = athena_read.athdf(file0)
    
    luminosity, B, time_list, entanglement = loop_func(fl, i_fl)

    luminosity = luminosity[1:]
    time_list  = time_list[1:]
    B = B[1:]
    entanglement = entanglement[1:]


    np.savetxt(f"save_arr/{save_arr_dir}luminosity_{fl}", np.array(luminosity))
    np.savetxt(f"save_arr/{save_arr_dir}time_list_{fl}", np.array(time_list))
    np.savetxt(f"save_arr/{save_arr_dir}B_{fl}", np.array(B))
    np.savetxt(f"save_arr/{save_arr_dir}entanglement_{fl}", np.array(entanglement))