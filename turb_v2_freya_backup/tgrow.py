import numpy as np


def tgrow (Chi, M, Rcl_lsh, Rcl_Lbox, tcool):

    t_grow =  Chi * np.sqrt(M) * np.sqrt(Rcl_lsh) * (Rcl_Lbox**(1/6)) * tcool

    return t_grow