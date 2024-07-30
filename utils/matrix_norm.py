from sympy import *
from sympy import Abs
import numpy as np
from numpy.linalg import eig


def norm_of_matrix(main_matrix):
    # To calculate the norm of a matrix M we use the Rayleigh-Ritz Theorem. Which assures us that
    # the operator norm or 2-norm of the matrix M is the square root of the largest eigenvalue of the matrix MT^{T}*M.

    mult_transpose_main = main_matrix.T * main_matrix
    eigenvalues_h, eigenvectors_h = eig(np.array(mult_transpose_main).astype(np.float64))

    # Now we calculate numerically the eigenvalues of M^{T}*M using numpy, which are non-negative.
    norm_eigenvalues_mtm = []
    for i in range(len(eigenvalues_h)):
        norm_eigenvalues_mtm.append(Abs(N(re(eigenvalues_h[i]), 10)))

    s = sqrt(max(norm_eigenvalues_mtm))
    return s
