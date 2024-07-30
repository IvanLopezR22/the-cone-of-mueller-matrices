import numpy as np
from sympy import *
from request_matrix import request_matrix
from appr_invertible_matrix import make_invertible
from appr_mueller_matrix import appr_mueller_matrix
from know_if_is_mueller import know_if_is_mueller

np.set_printoptions(precision=10, suppress=True, linewidth=2000)


def appr_invert_mueller_matrix(matrix_name: str, main_matrix: Matrix) -> Matrix:
    """
    Given a 4x4 matrix M, this function modify the diagonal of the matrix M so that the new matrix M(mu,inv,epsilon)
    is an invertible Mueller matrix.
    Using functions make_invertible(M), appr_mueller_matrix(M) and the notation of the comments from the files
    appr_invertible_matrix and appr_mueller_matrix we have 3 cases:
    Case 1: The matrix M is an invertible Mueller matrix, then M=M(mu,inv,epsilon).
    Case 2: The matrix M is a Mueller matrix but is not invertible, then M(mu,inv,epsilon)=M(inv,epsilon).
    Case 3: The matrix M is not a Mueller, then M(mu,inv,epsilon)=M(mu)(inv,epsilon).
    Parameters:
        -main_matrix (sympy 4x4 matrix): The matrix we want to approximate to an invertible Mueller matrix.

    Returns:
        Matrix: A sympy 4x4 matrix that is an invertible Mueller matrix.
    """
    m_mue_appr, name = appr_mueller_matrix(matrix_name, main_matrix)
    m_mue_inv, eps = make_invertible(name, m_mue_appr)
    qm, mi_qm, minimum_qm, first_column_is_stokes_m, fig = know_if_is_mueller(m_mue_inv)
    print(f"Due to the previous data, the approximation {matrix_name}(mu,inv,1/100) of {matrix_name} by an "
          f"invertible Mueller matrix is:\n{np.array(m_mue_inv).astype(np.float64)}")
    print(f"The minimum of the function q{matrix_name}(mu,inv,1/100) of {matrix_name}(mu,inv,1/100) is {minimum_qm}.")
    print(f"The first column of the matrix {matrix_name} is a Stokes vector: {first_column_is_stokes_m}.")
    print(f"The determinant of {matrix_name}(mu,inv,1/100) is {m_mue_inv.det()}.\nTherefore, "
          f"{matrix_name}(mu,inv,1/100) is an invertible Mueller matrix.")

    return m_mue_inv


if __name__ == '__main__':
    m_0 = request_matrix("M")
    print("------------------------------------------------------------------")
    print(f"The input matrix M is:\n{np.array(m_0).astype(np.float64)}")
    print("------------------------------------------------------------------")
    appr_invert_mueller_matrix("M", m_0)
