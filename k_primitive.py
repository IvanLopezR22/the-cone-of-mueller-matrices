from sympy import *
import numpy as np
from numpy.linalg import eig
from know_if_is_mueller import know_if_is_mueller
from utils.request_matrix import request_matrix


def k_primitive(matrix_name: str, main_matrix: Matrix):
    """
              This functions uses the fact that a Matrix M is K-irreducible if and only if:
              1. The minimum of the function qM of M is greater than zero;.
              2. The first column of the matrix M is a Stokes vector.
               Parameters:
              -matrix_name (String): The name of the matrix we are working with.
              -main_matrix (sympy 4x4 matrix): The matrix to which we want to calculate the norm.

              Returns:
              Nothing.
          """

    # The matrix is K-primitive if and only if the minimums of the function qM is positive and the first
    # column of the matrix is a Stokes vector.
    qm, mi_qm, minimum_qm, first_column_is_stokes, fig = know_if_is_mueller(main_matrix)
    if minimum_qm > 0 and first_column_is_stokes == True:
        print(f"The minimum of the function q{matrix_name} of the matrix {matrix_name} is {minimum_qm}.")
        print(f"The first column of the matrix {matrix_name} is a Stokes vector: {first_column_is_stokes}.")
        print(f"Therefore, the matrix {matrix_name} is a K-primitive Mueller matrix.")
    else:
        print(f"The minimum of the function q{matrix_name} of the matrix {matrix_name} is {minimum_qm}.")
        print(f"The first column of the matrix {matrix_name} is a Stokes vector: {first_column_is_stokes}.")
        print(f"Therefore, the matrix {matrix_name} is NOT a K-primitive Mueller matrix")



