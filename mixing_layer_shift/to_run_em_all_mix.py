import numpy as np
import os

import para_scan as ps 
import create_input as ci
import create_job_script as cj


for i in range(np.size(ps.box_width)):
    for j in range(np.size(ps.Ma)):
        for k in range(np.size(ps.B_dir)):

             ## os.system('source ../load_module')
            os.system(f'cp -r template_dir mix'+ps.filename_mix_add(i,j,k))

            os.system(f'rm mix'+ps.filename_mix_add(i,j,k)+'/job_script_mix_template.sh')
            os.system(f'rm mix'+ps.filename_mix_add(i,j,k)+'/athinput_mix_template')

            print(f"\n\n{ps.filename_mix_add(i,j,k)}: Directory created...")
            ci.input_mix_new(i,j,k)
            print("Input file created....")
            cj.jobscript_mix_new(i,j,k)
            print("Job script created....")

            if ps.B_flag:
                os.system(f'(cd mix{ps.filename_mix_add(i,j,k)} && ./turb_compile_freya_wB)')
                # os.system(f'(cd para_scan{filename_add(i,j)} && ./turb_compile_wB)')
            else:
                os.system(f'(cd mix{ps.filename_mix_add(i,j,k)} && ./turb_compile_freya)')
                # os.system(f'(cd para_scan{filename_add(i,j)} && ./turb_compile)')

            # os.system(f'(cd mix{ps.filename_mix_add(i,j,k)} && sbatch job_script_mix{ps.filename_mix_add(i,j,k)}.sh)')
            # os.system(f'(cd para_scan{filename_add(i,j)} && ./turb_run athinput{filename_add(i,j)}.turb 2)')

            print("_________________________________________________")
