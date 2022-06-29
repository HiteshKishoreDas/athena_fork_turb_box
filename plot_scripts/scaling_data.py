import numpy as np
from scaling_class import scaling_data

scale_neg_rseed = []
scale_pos_rseed = []


# Negative rseed runs

## 64^3 boxes

# scale_neg_rseed.append(scaling_data(axis_cells=[64,64,64], dim=3, nproc=64, time_sec=11*3600+25*60+33, \
#     mesh_block=[16,16,16]))


## 128^3 boxes

scale_neg_rseed.append(scaling_data(axis_cells=[128,128,128], dim=3, nproc=64, time_sec=6*60+39, \
    mesh_block=[32,32,32]))

## 256^3 boxes

scale_neg_rseed.append(scaling_data(axis_cells=[256,256,256], dim=3, nproc=64*8, time_sec=15*60+26, \
    mesh_block=[32,32,32]))


# Postive rseed runs

## 64^3 boxes

# scale_pos_rseed.append(scaling_data(axis_cells=[64,64,64], dim=3, nproc=64, time_sec=11*3600+21*60+25, \
#     mesh_block=[16,16,16],pos_rseed=True))


## 128^3 boxes


if __name__=="__main__":

    if len(scale_neg_rseed)!=0:
        scale_neg_rseed[len(scale_neg_rseed)-1].show_data()
        
    print('\n')

    if len(scale_pos_rseed)!=0:
        scale_pos_rseed[len(scale_pos_rseed)-1].show_data()