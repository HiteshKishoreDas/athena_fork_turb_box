import numpy as np
import para_scan as ps

def tag_replace(filedata,tag,value):

    filedata = filedata.replace(tag,str(value))

    return filedata

def athinput_replace(template_data,i,j,k):

    template_data = tag_replace(template_data,"!!HST_DT!!" ,ps.hst_dt_arr[i])
    template_data = tag_replace(template_data,"!!OUT_DT!!" ,ps.out_dt_arr[i])
    template_data = tag_replace(template_data,"!!RST_DT!!" ,ps.rst_dt_arr[i])
    template_data = tag_replace(template_data,"!!USER_OUT_DT!!",ps.out_dt_arr[i])
    template_data = tag_replace(template_data,"!!TLIM!!"   ,ps.tlim[i])


    template_data = tag_replace(template_data,"!!NX1!!",ps.nx1[0])
    template_data = tag_replace(template_data,"!!NX2!!",ps.nx2[0])
    template_data = tag_replace(template_data,"!!NX3!!",ps.nx3[0])

    template_data = tag_replace(template_data,"!!NX1_MESH!!",ps.nx1_mesh[0])
    template_data = tag_replace(template_data,"!!NX2_MESH!!",ps.nx2_mesh[0])
    template_data = tag_replace(template_data,"!!NX3_MESH!!",ps.nx3_mesh[0])

    template_data = tag_replace(template_data,"!!X1MIN!!",ps.x1min[i])
    template_data = tag_replace(template_data,"!!X2MIN!!",ps.x2min[i])
    template_data = tag_replace(template_data,"!!X3MIN!!",ps.x3min[i])

    template_data = tag_replace(template_data,"!!X1MAX!!",ps.x1max[i])
    template_data = tag_replace(template_data,"!!X2MAX!!",ps.x2max[i])
    template_data = tag_replace(template_data,"!!X3MAX!!",ps.x3max[i])

    template_data = tag_replace(template_data,"!!P_FLOOR!!",ps.P_floor[0])

    template_data = tag_replace(template_data,"!!AMB_RHO!!",ps.amb_rho[0])

    template_data = tag_replace(template_data,"!!FRONT_TH!!",ps.front_thickness[i])
    template_data = tag_replace(template_data,"!!V_SHEAR!!",ps.v_shear)
    template_data = tag_replace(template_data,"!!V_SHIFT!!",ps.v_shift)

    template_data = tag_replace(template_data,"!!KNX_KH!!",ps.knx_KH)
    template_data = tag_replace(template_data,"!!KNY_KH!!",ps.kny_KH)
    template_data = tag_replace(template_data,"!!AMP_KH!!",ps.amp_KH)

    template_data = tag_replace(template_data,"!!T_FLOOR!!"  ,ps.T_floor)
    template_data = tag_replace(template_data,"!!T_CEIL!!"   ,ps.T_ceil)
    template_data = tag_replace(template_data,"!!T_HOT!!"    ,ps.T_hot)
    template_data = tag_replace(template_data,"!!T_COLD!!"   ,ps.T_cold)
    template_data = tag_replace(template_data,"!!T_CUT_MUL!!",ps.T_cut_mul)

    template_data = tag_replace(template_data,"!!XSOL!!",ps.Xsol)
    template_data = tag_replace(template_data,"!!ZSOL!!",ps.Zsol)

    template_data = tag_replace(template_data,"!!LAM_FAC!!",ps.Lambda_fac[0])

    template_data = tag_replace(template_data,"!!B_X!!",ps.B_x[i,j,k])
    template_data = tag_replace(template_data,"!!B_Y!!",ps.B_y[i,j,k])
    template_data = tag_replace(template_data,"!!B_Z!!",ps.B_z[i,j,k])


    return template_data


def input_mix_new(i, j, k, template='template_dir/athinput_mix_template'):

    input_newfile = f'mix'+ps.filename_mix_add(i,j,k)+'/athinput_mix'+ps.filename_mix_add(i,j,k)

    with open(template,'r') as file:
        template_data = file.read()

    template_data = athinput_replace(template_data,i,j,k)

    with open(input_newfile,'w') as file:
        file.write(template_data)
