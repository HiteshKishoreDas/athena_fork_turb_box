import numpy as np
import os

from turb_template.template_dir.para_scan import *
import turb_template.template_dir.create_input as ci
import turb_template.template_dir.create_job_script as cj


# for i in range(np.size(R_lsh)): # np.size(nx1)):
# for i in range(np.size(R_lsh)): # np.size(nx1)):
for i in [len(R_lsh) - 1]:  # np.size(nx1)):
    for j in range(np.size(nx1)):
        ## os.system('source ../load_module')
        os.system(f"cp -r template_dir para_scan" + filename_turb_add(i, j))

        os.system(
            f"rm para_scan" + filename_turb_add(i, j) + "/job_script_turb_template.sh"
        )
        os.system(
            f"rm para_scan"
            + filename_turb_add(i, j)
            + "/athinput.turb_v2_turb_template"
        )

        os.system(
            f"rm para_scan" + filename_turb_add(i, j) + "/job_script_cloud_template.sh"
        )
        os.system(
            f"rm para_scan"
            + filename_turb_add(i, j)
            + "/athinput.turb_v2_cloud_template"
        )

        print(f"\n\n{filename_turb_add(i,j)}: Directory created...")
        ci.input_turb_new(i, j)
        print("Input file created....")
        cj.jobscript_turb_new(i, j)
        print("Job script created....")

        if B_flag:
            os.system(
                f"(cd para_scan{filename_turb_add(i,j)} && ./turb_compile_freya_wB)"
            )
            # os.system(f'(cd para_scan{filename_add(i,j)} && ./turb_compile_wB)')
        else:
            os.system(f"(cd para_scan{filename_turb_add(i,j)} && ./turb_compile_freya)")
            # os.system(f'(cd para_scan{filename_add(i,j)} && ./turb_compile)')

        os.system(
            f"(cd para_scan{filename_turb_add(i,j)} && sbatch job_script_turb{filename_turb_add(i,j)}.sh)"
        )
        # os.system(f'(cd para_scan{filename_add(i,j)} && ./turb_run athinput{filename_add(i,j)}.turb 2)')

        print("_________________________________________________")
