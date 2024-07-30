from sympy import *
from sympy import Abs
import numpy as np
from numpy.linalg import eig
from request_matrix import request_matrix

np.set_printoptions(precision=10, suppress=True, linewidth=2000)


def matrix_norm(main_matrix: Matrix) -> Matrix:
    """
        This functions calculates the operator norm or 2-norm of a 4x4 matrix with real entries.
        The 2-norm is defined as ||M||_2=sup{||M(x)|| : ||x||=1}.
        To calculate the norm of a matrix M we use the Rayleigh-Ritz Theorem. Which assures us that
        the operator norm or 2-norm of the matrix M is the square root of the largest eigenvalue of the matrix MT^{T}*M.

         Parameters:
        -main_matrix (sympy 4x4 matrix): The matrix to which we want to calculate the norm.

        Returns:
        Float: The norm of the matrix, a positive real number.
    """

    mult_transpose_main = main_matrix.T * main_matrix
    eigenvalues_h, eigenvectors_h = eig(np.array(mult_transpose_main).astype(np.float64))

    # Now we calculate numerically the eigenvalues of M^{T}*M using numpy, which are non-negative.
    norm_eigenvalues_mtm = []
    for i in range(len(eigenvalues_h)):
        norm_eigenvalues_mtm.append(Abs(N(re(eigenvalues_h[i]), 10)))

    s = sqrt(max(norm_eigenvalues_mtm))
    return s


if __name__ == '__main__':
    m = request_matrix("M")
    print(f"The input matrix M is:\n{np.array(m).astype(np.float64)}")
    print("----------------------------------------")
    print(f"The norm of the introduced matrix M is: {matrix_norm(m)}")
    print("----------------------------------------")
