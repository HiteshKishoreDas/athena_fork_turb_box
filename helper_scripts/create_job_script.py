import numpy as np
import turb_template.template_dir.para_scan as ps


def tag_replace(filedata, tag, value):
    filedata = filedata.replace(tag, str(value))

    return filedata


def jobscript_turb_new(i, j, template="template_dir/job_script_turb_template.sh"):
    jobscript_newfile = (
        f"para_scan"
        + ps.filename_turb_add(i, j)
        + "/job_script_turb"
        + ps.filename_turb_add(i, j)
        + ".sh"
    )

    with open(template, "r") as file:
        template_data = file.read()

    if ps.B_flag:
        # job_name = f'm{i}{j}{ps.rseed}'
        job_name = f"m{i}{ps.rseed}{ps.M}"
    else:
        # job_name = f'h{i}{j}{ps.rseed}'
        job_name = f"h{i}{ps.rseed}{ps.M}"

    work_dir = ps.parent_dir + f"/para_scan" + ps.filename_turb_add(i, j) + "/"
    input_file = f"athinput_turb" + ps.filename_turb_add(i, j) + ".turb"

    template_data = tag_replace(template_data, "!!JOB_NAME!!", job_name)
    template_data = tag_replace(template_data, "!!QUEUE!!", ps.queue)
    template_data = tag_replace(template_data, "!!NODES!!", ps.nodes[j])

    template_data = tag_replace(
        template_data, "!!NTASKS_PER_NODE!!", ps.ntasks_per_node
    )
    template_data = tag_replace(template_data, "!!TIME_LIMIT!!", ps.time_limit_turb[j])
    template_data = tag_replace(
        template_data, "!!TIME_LIMIT_RST!!", ps.time_limit_turb_rst[j]
    )

    template_data = tag_replace(template_data, "!!DIR_PATH_ADD!!", ps.dir_path_add)

    template_data = tag_replace(template_data, "!!WORK_DIR!!", work_dir)
    template_data = tag_replace(template_data, "!!INPUT_FILE!!", input_file)

    with open(jobscript_newfile, "w") as file:
        file.write(template_data)


def jobscript_cloud_new(i, j, template="template_dir/job_script_cloud_template.sh"):
    jobscript_newfile = (
        f"para_scan"
        + ps.filename_cloud_add(i, j)
        + "/job_script_cloud"
        + ps.filename_cloud_add(i, j)
        + ".sh"
    )

    with open(template, "r") as file:
        template_data = file.read()

    if ps.B_flag:
        job_name = f"m{i}{j}M{ps.i_mach}{ps.rseed}"
    else:
        job_name = f"h{i}{j}M{ps.i_mach}{ps.rseed}"

    work_dir = f"{ps.parent_dir}/para_scan{ps.filename_cloud_add(i,j)}/"
    input_file = f"athinput_cloud{ps.filename_cloud_add(i,j)}.turb"
    turb_dir = f"para_scan{ps.filename_turb_add(i,j)}"
    rst_file = f"../Turbulence/{turb_dir}/Turb.final.rst"

    template_data = tag_replace(template_data, "!!JOB_NAME!!", job_name)
    template_data = tag_replace(template_data, "!!QUEUE!!", ps.queue)
    template_data = tag_replace(template_data, "!!NODES!!", ps.nodes[j])

    template_data = tag_replace(
        template_data, "!!NTASKS_PER_NODE!!", ps.ntasks_per_node
    )
    template_data = tag_replace(template_data, "!!TIME_LIMIT!!", ps.time_limit_cloud[j])
    template_data = tag_replace(
        template_data, "!!TIME_LIMIT_RST!!", ps.time_limit_cloud_rst[j]
    )

    template_data = tag_replace(template_data, "!!WORK_DIR!!", work_dir)
    template_data = tag_replace(template_data, "!!TURB_DIR!!", turb_dir)

    template_data = tag_replace(template_data, "!!DIR_PATH_ADD!!", ps.dir_path_add)

    template_data = tag_replace(template_data, "!!INPUT_FILE!!", input_file)
    template_data = tag_replace(template_data, "!!RST_FILE!!", rst_file)

    with open(jobscript_newfile, "w") as file:
        file.write(template_data)
