import numpy as np
import os

import para_scan as ps
import create_input as ci
import create_job_script as cj


# for i in range(np.size(R_lsh)): # np.size(nx1)):
# for i in range(np.size(R_lsh)): # np.size(nx1)):
# for i in [len(R_lsh)-1]: # np.size(nx1)):
for i in [1, 2]:
    for j in range(np.size(ps.nx1)):

        # os.system("module purge")
        # os.system("source ../load_module_cobra")
        os.system(
            f"cp -r template_dir_turb Turbulence/para_scan" + ps.filename_turb_add(i, j)
        )
        print("Template directory copied...")

        os.system(
            f"rm Turbulence/para_scan"
            + ps.filename_turb_add(i, j)
            + "/job_script_turb_template.sh"
        )
        print("Template job script removed...")

        os.system(
            f"rm Turbulence/para_scan"
            + ps.filename_turb_add(i, j)
            + "/athinput.turb_v2_turb_template"
        )
        print("Template input file removed...")

        print(f"\n\n{ps.filename_turb_add(i,j)}: Directory created...")
        ci.input_turb_new(i, j)
        print("Input file created....")
        cj.jobscript_turb_new(i, j)
        print("Job script created....")

        if ps.B_flag:
            os.system(
                f"(cd Turbulence/para_scan{ps.filename_turb_add(i,j)} && ./turb_compile_freya_wB)"
            )
            # os.system(f'(cd para_scan{filename_add(i,j)} && ./turb_compile_wB)')
        else:
            os.system(
                f"(cd Turbulence/para_scan{ps.filename_turb_add(i,j)} && ./turb_compile_freya)"
            )
            # os.system(f'(cd para_scan{filename_add(i,j)} && ./turb_compile)')

        # os.system(
        #     f"(cd para_scan{filename_turb_add(i,j)} && sbatch job_script_turb{filename_turb_add(i,j)}.sh)"
        # )
        # os.system(f'(cd para_scan{filename_add(i,j)} && ./turb_run athinput{filename_add(i,j)}.turb 2)')

        print("_________________________________________________")
