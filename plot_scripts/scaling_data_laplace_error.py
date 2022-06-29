import numpy as np
from scaling_class import scaling_data

scale_neg_rseed = []
scale_pos_rseed = []


# Negative rseed runs

## 32^3 boxes

scale_neg_rseed.append(scaling_data(axis_cells=[32,32,32], dim=3, nproc=1, time_sec=199.451))

scale_neg_rseed.append(scaling_data(axis_cells=[32,32,32], dim=3, nproc=2, time_sec=133.524, \
    mesh_block=[16,32,32]))


## 64^3 boxes

scale_neg_rseed.append(scaling_data(axis_cells=[64,64,64], dim=3, nproc=8, time_sec=5*60+35, \
    mesh_block=[32,32,32]))

scale_neg_rseed.append(scaling_data(axis_cells=[64,64,64], dim=3, nproc=16, time_sec=4*60+12, \
    mesh_block=[16,32,32]))

scale_neg_rseed.append(scaling_data(axis_cells=[64,64,64], dim=3, nproc=32, time_sec=3*60+25, \
    mesh_block=[16,16,32]))

scale_neg_rseed.append(scaling_data(axis_cells=[64,64,64], dim=3, nproc=64, time_sec=1*60+53, \
    mesh_block=[16,16,16]))


## 128^3 boxes

scale_neg_rseed.append(scaling_data(axis_cells=[128,128,128], dim=3, nproc=8, time_sec=1*60*60+15*60+7, \
    mesh_block=[64,64,64]))

scale_neg_rseed.append(scaling_data(axis_cells=[128,128,128], dim=3, nproc=16, time_sec=48*60+7, \
    mesh_block=[32,64,64]))

scale_neg_rseed.append(scaling_data(axis_cells=[128,128,128], dim=3, nproc=32, time_sec=25*60+58, \
    mesh_block=[32,32,64]))

scale_neg_rseed.append(scaling_data(axis_cells=[128,128,128], dim=3, nproc=64, time_sec=15*60+43, \
    mesh_block=[32,32,32]))



if __name__=="__main__":

    scale_neg_rseed[len(scale_neg_rseed)-1].show_data()