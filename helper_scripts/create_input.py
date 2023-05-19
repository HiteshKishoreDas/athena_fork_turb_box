import numpy as np
import turb_template.template_dir.para_scan as ps


def tag_replace(filedata, tag, value):
    filedata = filedata.replace(tag, str(value))

    return filedata


def athinput_replace(template_data, i, j):
    template_data = tag_replace(template_data, "!!HST_DT_TURB!!", ps.hst_dt_trb_arr[i])
    template_data = tag_replace(template_data, "!!OUT_DT_TURB!!", ps.out_dt_trb_arr[i])
    template_data = tag_replace(template_data, "!!RST_DT_TURB!!", ps.rst_dt_trb_arr[i])
    template_data = tag_replace(template_data, "!!USER_DT_TURB!!", ps.out_dt_trb_arr[i])
    template_data = tag_replace(template_data, "!!TLIM_TURB!!", ps.tlim_trb[i])

    template_data = tag_replace(template_data, "!!HST_DT_CLD!! ", ps.hst_dt_cld_arr[i])
    template_data = tag_replace(template_data, "!!OUT_DT_CLD!! ", ps.out_dt_cld_arr[i])
    template_data = tag_replace(template_data, "!!RST_DT_CLD!! ", ps.rst_dt_cld_arr[i])
    template_data = tag_replace(template_data, "!!USER_DT_CLD!!", ps.out_dt_cld_arr[i])
    template_data = tag_replace(template_data, "!!TLIM_CLD!!", ps.tlim_cld[i])

    template_data = tag_replace(template_data, "!!NX1!!", ps.nx1[j])
    template_data = tag_replace(template_data, "!!NX2!!", ps.nx2[j])
    template_data = tag_replace(template_data, "!!NX3!!", ps.nx3[j])

    template_data = tag_replace(template_data, "!!NX1_MESH!!", ps.nx1_mesh[j])
    template_data = tag_replace(template_data, "!!NX2_MESH!!", ps.nx2_mesh[j])
    template_data = tag_replace(template_data, "!!NX3_MESH!!", ps.nx3_mesh[j])

    template_data = tag_replace(template_data, "!!X1MIN!!", ps.x1min[i])
    template_data = tag_replace(template_data, "!!X2MIN!!", ps.x2min[i])
    template_data = tag_replace(template_data, "!!X3MIN!!", ps.x3min[i])

    template_data = tag_replace(template_data, "!!X1MAX!!", ps.x1max[i])
    template_data = tag_replace(template_data, "!!X2MAX!!", ps.x2max[i])
    template_data = tag_replace(template_data, "!!X3MAX!!", ps.x3max[i])

    template_data = tag_replace(template_data, "!!P_FLOOR!!", ps.P_floor)

    template_data = tag_replace(template_data, "!!DEDT!!", ps.dedt[i])

    template_data = tag_replace(template_data, "!!T_CORR!!", ps.t_corr_arr[i])
    template_data = tag_replace(template_data, "!!DT_DRIVE!!", ps.dt_drive_arr[i])
    template_data = tag_replace(template_data, "!!F_SHEAR!!", ps.f_shear)
    template_data = tag_replace(template_data, "!!RSEED!!", ps.rseed)

    # template_data = tag_replace(template_data,"!!COOL_FLAG!!" ,ps.cooling_flag)
    # template_data = tag_replace(template_data,"!!CLOUD_FLAG!!",ps.cloud_flag)

    template_data = tag_replace(template_data, "!!AMB_RHO!!", ps.amb_rho)

    template_data = tag_replace(template_data, "!!CLOUD_RADIUS!!", ps.cloud_radius[i])
    template_data = tag_replace(template_data, "!!CLOUD_TIME!!", ps.cloud_time[i])
    template_data = tag_replace(template_data, "!!CLOUD_CHI!!", ps.cloud_chi)

    template_data = tag_replace(template_data, "!!T_FLOOR!!", ps.T_floor)
    template_data = tag_replace(template_data, "!!T_CEIL!!", ps.T_ceil)
    template_data = tag_replace(template_data, "!!T_HOT!!", ps.T_hot)
    template_data = tag_replace(template_data, "!!T_HOT_REQ!!", ps.T_hot_req)
    template_data = tag_replace(template_data, "!!T_COLD!!", ps.T_cold)
    template_data = tag_replace(template_data, "!!T_CUT_MUL!!", ps.T_cut_mul)
    template_data = tag_replace(template_data, "!!T_CUT!!", ps.T_cut)

    template_data = tag_replace(template_data, "!!XSOL!!", ps.Xsol)
    template_data = tag_replace(template_data, "!!ZSOL!!", ps.Zsol)

    template_data = tag_replace(template_data, "!!CLOUD_POS_X!!", ps.cloud_pos_x)
    template_data = tag_replace(template_data, "!!CLOUD_POS_Y!!", ps.cloud_pos_y)
    template_data = tag_replace(template_data, "!!CLOUD_POS_Z!!", ps.cloud_pos_z)

    if ps.B_flag:
        template_data = tag_replace(template_data, "!!B_X!!", ps.B_x)
        template_data = tag_replace(template_data, "!!B_Y!!", ps.B_y)
        template_data = tag_replace(template_data, "!!B_Z!!", ps.B_z)

    else:
        template_data = tag_replace(template_data, "!!B_X!!", 0.0)
        template_data = tag_replace(template_data, "!!B_Y!!", 0.0)
        template_data = tag_replace(template_data, "!!B_Z!!", 0.0)

    return template_data


def input_turb_new(i, j, template="template_dir/athinput.turb_v2_turb_template"):
    input_newfile = (
        f"para_scan"
        + ps.filename_turb_add(i, j)
        + "/athinput_turb"
        + ps.filename_turb_add(i, j)
        + ".turb"
    )

    with open(template, "r") as file:
        template_data = file.read()

    template_data = athinput_replace(template_data, i, j)

    with open(input_newfile, "w") as file:
        file.write(template_data)


def input_cloud_new(i, j, template="template_dir/athinput.turb_v2_cloud_template"):
    input_newfile = (
        f"para_scan"
        + ps.filename_cloud_add(i, j)
        + "/athinput_cloud"
        + ps.filename_cloud_add(i, j)
        + ".turb"
    )

    with open(template, "r") as file:
        template_data = file.read()

    template_data = athinput_replace(template_data, i, j)

    with open(input_newfile, "w") as file:
        file.write(template_data)
