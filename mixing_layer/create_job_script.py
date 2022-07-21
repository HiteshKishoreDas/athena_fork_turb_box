import numpy as np
import para_scan as ps 

def tag_replace(filedata,tag,value):

    filedata = filedata.replace(tag,str(value))

    return filedata


def jobscript_mix_new(i, j, k, template='template_dir/job_script_mix_template.sh'):

    jobscript_newfile = f'mix'+ps.filename_mix_add(i,j,k)+'/job_script_mix'+ps.filename_mix_add(i,j,k)+'.sh'

    with open(template,'r') as file:
        template_data = file.read()

    if ps.B_flag:
        job_name = f'm{i}{j}{k}'
    else:
        job_name = f'h{i}{j}{k}'
    
    work_dir = ps.parent_dir+f'/mix'+ps.filename_mix_add(i,j,k)+'/'
    input_file = f'athinput_mix'+ps.filename_mix_add(i,j,k)+'.turb'

    template_data = tag_replace(template_data,"!!JOB_NAME!!",job_name)
    template_data = tag_replace(template_data,"!!QUEUE!!",ps.queue)
    template_data = tag_replace(template_data,"!!NODES!!",ps.nodes[0])

    template_data = tag_replace(template_data,"!!NTASKS_PER_NODE!!",ps.ntasks_per_node)
    template_data = tag_replace(template_data,"!!TIME_LIMIT!!",ps.time_limit[0])
    template_data = tag_replace(template_data,"!!TIME_LIMIT_RST!!",ps.time_limit_rst[0])


    template_data = tag_replace(template_data,"!!WORK_DIR!!",work_dir)
    template_data = tag_replace(template_data,"!!INPUT_FILE!!",input_file)

    with open(jobscript_newfile,'w') as file:
        file.write(template_data)