from sympy import *
import numpy as np
from numpy.linalg import eig
from know_if_is_mueller import know_if_is_mueller
from utils.request_matrix import request_matrix


def k_primitive(m):
    # The matrix is K-primitive if and only if the minimums of the functions qM and proj1(M) are positive.
    qm, mi_qm, minimum_qm, proy1_m, mi_proy1_m, minimum_proy1_m, fig = know_if_is_mueller(m)
    if minimum_qm > 0 and minimum_proy1_m > 0:
        return True
    else:
        return False
    