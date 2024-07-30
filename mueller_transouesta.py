from appr_mueller_matrix import appr_mueller_matrix
from utils.request_matrix import request_matrix_interface
from know_if_is_mueller import know_if_is_mueller
from sympy import *

m = request_matrix_interface("M")
m_appr, name, is_mueller = appr_mueller_matrix("M", m)
print(f"Our matrix M is:")
print(m_appr)
#m_appr = Matrix([[83.4711680262525, 4, 6, 5], [2, 1, 5, 2], [3, 6, 9, 85], [6, 3, 2, 1]])
qm, mi_qm, minimum_qm, proy1_m, mi_proy1_m, minimum_proy1_m, fig = know_if_is_mueller(m_appr)

if(minimum_qm>=0 and minimum_proy1_m>=0):
    print(f"The minimum of the functions qM and proy1(M) of the matrix M are "
          f"{minimum_qm} and {minimum_proy1_m} respectively.")
    print(f"The matrix M is Mueller")
else:
    print(f"The minimum of the functions qM and proy1(M) of the matrix M are "
          f"{minimum_qm} and {minimum_proy1_m} respectively.")
    print(f"The matrix M is NOT Mueller")

qm_1, mi_qm_1, minimum_qm_1, proy1_m_1, mi_proy1_m_1, minimum_proy1_m_1, fig_1 = know_if_is_mueller(m_appr.T)

if(minimum_qm_1>=0 and minimum_proy1_m_1>=0):
    print(f"The minimum of the functions qM and proy1(M) of the matrix M are "
          f"{minimum_qm_1} and {minimum_proy1_m_1} respectively.")
    print(f"The transpose of the matrix M is Mueller")
else:
    print(f"The minimum of the functions qM and proy1(M) of the matrix M are "
          f"{minimum_qm_1} and {minimum_proy1_m_1} respectively.")
    print(f"The transpose of the matrix M is NOT Mueller")




