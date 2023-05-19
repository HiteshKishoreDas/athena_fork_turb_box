import numpy as np
import os

from turb_template.template_dir.para_scan import *
import turb_template.template_dir.create_input as ci
import turb_template.template_dir.create_job_script as cj


def job_list_read(filename):
    job_list = []

    # reading csv file
    with open(filename, "r") as fl:
        end = False

        while not (end):
            l = fl.readline()
            if not l:
                end = True
                break
            # print(l)
            job_list.append(
                [
                    l.split(" ")[0],
                    l.split(" ")[1][0:1],
                    l.split(" ")[1][1:2],
                    0,
                    l.split(" ")[1][2:3],
                    l.split(" ")[1][3:],
                ]
            )

    return job_list


def job_id(job_list, i, j, B_flag):
    if B_flag:
        # mhd = 'mhd'
        mhd = "m"
    else:
        # mhd = 'hyd'
        mhd = "h"

    for job in job_list:
        if job[1] == mhd:
            if i == int(job[2]):
                if j == int(job[3]):
                    return job[0]

    return -1


job_list = job_list_read("job_list/M_0.5_job_list")


# for i in [4]:#range(np.size(R_lsh)):
# for i in range(np.size(R_lsh)):
for i in [len(R_lsh) - 1]:  # range(np.size(R_lsh)):
    for j in range(np.size(nx1)):
        turb_job_id = job_id(job_list, i, j, B_flag)

        if turb_job_id == -1:
            continue

        ## os.system('source ../load_module')
        os.system(f"cp -r template_dir para_scan" + filename_cloud_add(i, j))

        os.system(
            f"rm para_scan" + filename_cloud_add(i, j) + "/job_script_turb_template.sh"
        )
        os.system(
            f"rm para_scan"
            + filename_cloud_add(i, j)
            + "/athinput.turb_v2_turb_template"
        )

        os.system(
            f"rm para_scan" + filename_cloud_add(i, j) + "/job_script_cloud_template.sh"
        )
        os.system(
            f"rm para_scan"
            + filename_cloud_add(i, j)
            + "/athinput.turb_v2_cloud_template"
        )

        # os.system(f'cp para_scan{filename_turb_add(i,j)}/Turb.hst para_scan{filename_cloud_add(i,j)}/')

        print(f"\n\n{filename_cloud_add(i,j)}: Directory created...")
        ci.input_cloud_new(i, j)
        print("Input file created....")
        cj.jobscript_cloud_new(i, j)
        print("Job script created....")

        if B_flag:
            os.system(
                f"(cd para_scan{filename_cloud_add(i,j)} && ./turb_compile_freya_wB)"
            )
            # os.system(f'(cd para_scan{filename_add(i,j)} && ./turb_compile_wB)')
        else:
            os.system(
                f"(cd para_scan{filename_cloud_add(i,j)} && ./turb_compile_freya)"
            )
            # os.system(f'(cd para_scan{filename_add(i,j)} && ./turb_compile)')

        os.system(
            f"(cd para_scan{filename_cloud_add(i,j)} && sbatch --dependency=afterok:{turb_job_id} job_script_cloud{filename_cloud_add(i,j)}.sh)"
        )
        # os.system(f'(cd para_scan{filename_cloud_add(i,j)} && sbatch job_script_cloud{filename_cloud_add(i,j)}.sh)')

        # os.system(f'(cd para_scan{filename_add(i,j)} && ./turb_run athinput{filename_add(i,j)}.turb 2)')

        print("_________________________________________________")
