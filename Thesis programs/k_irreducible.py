from sympy import *
import numpy as np
from numpy.linalg import eig
from know_if_is_mueller import know_if_is_mueller
from request_matrix import request_matrix

np.set_printoptions(precision=10, suppress=True, linewidth=2000)


def k_irreducible(matrix_name: str, main_matrix: Matrix):
    """
           This functions uses the fact that a Matrix M is K-irreducible if and only if:
           1. The spectral radius p(M) of M is a simple eigenvalue.
           2. Let v the eigenvector associated to p(M). Then v or -v is in the interior of the Stokes cone.
           3. Any other eigenvector of M is not a Stokes vector.

            Parameters:
           -matrix_name (String): The name of the matrix we are working with.
           -main_matrix (sympy 4x4 matrix): The matrix to which we want to calculate the norm.

           Returns:
           Nothing.
       """

    # First we need to know if the matrix is Mueller.
    qm, mi_qm, minimum_qm, first_column_is_stokes_m, fig = know_if_is_mueller(main_matrix)
    if minimum_qm >= 0 and first_column_is_stokes_m == True:
        print(f"The minimum of the function q{matrix_name} of the matrix {matrix_name} is {minimum_qm}.")
        print(f"The first column of the matrix {matrix_name} is a Stokes vector: {first_column_is_stokes_m}.")
        print(f"Therefore the matrix {matrix_name} is a Mueller matrix..")
    else:
        print(f"The minimum of the function q{matrix_name} of the matrix {matrix_name} is {minimum_qm}.")
        print(f"The first column of the matrix {matrix_name} is a Stokes vector: {first_column_is_stokes_m}.")
        print(f"Therefore the matrix {matrix_name} is NOT a Mueller matrix. Hence, M is NOT K-irreducible.")
        return

    m_0 = np.array(main_matrix).astype(np.float64)
    # "eigenvalues" is the list of eigenvalues of the matrix M.
    # "eigenvectors" is a matrix that contains the eigenvectors in the columns.
    # The list of eigenvalues and the matrix of eigenvectors are related by the index, that is, the element i of the
    # list of eigenvalues corresponds the column i of the matrix of eigenvectors.
    eigenvalues, eigenvectors = eig(m_0)

    # We have to calculate the spectral radius of M.

    # When the eigenvalues are approximated numerically, the numerical approximations generate small imaginary parts
    # in the real eigenvalues. Then we filter that error into the following list, using the fact that in the case of
    # polynomials with real coefficients, if a complex number is a root, so is its conjugate.
    eigenvalues_without_error = []
    for i in eigenvalues:
        conj = i.conjugate()
        z = conj in eigenvalues
        if simplify(i).is_real == True:
            eigenvalues_without_error.append(np.round(i, 12))
        elif simplify(i).is_real == False and z == False:
            eigenvalues_without_error.append(np.round(i.real, 12))
        else:
            eigenvalues_without_error.append(i)
    print(f"The eigenvalues of the matrix {matrix_name} are: {eigenvalues_without_error}.")
    print(f"The eigenvectors (ordered in the columns respectively) of the matrix {matrix_name} are:\n{eigenvectors}")

    # norm_eigenvalues is the list of the norms of the eigenvalues.
    norm_eigenvalues = (np.absolute(eigenvalues_without_error)).tolist()
    print(f"The norms of the eigenvalues of {matrix_name} are: {norm_eigenvalues}.")

    # The spectral radius ro_m is the maximum in the list norm_eigenvalue.
    rho_m = max(norm_eigenvalues)

    # Now we ask if ro_m is in the list of eigenvalues_without_error.
    if rho_m in eigenvalues_without_error:
        print(f"The spectral radius, {rho_m}, is an eigenvalue of M.")
        ind_rho_m = eigenvalues_without_error.index(rho_m)
    else:
        print(f"The spectral radius, {rho_m}, is not an eigenvalue of M.")
        print(f"Therefore, the matrix {matrix_name} is NOT K-irreducible.")
        return

    # We need to know if the spectral radius is a simple eigenvalue. That is, we need to know if the algebraic
    # multiplicity of the spectral radius is 1.
    mult_alg_rho_m = eigenvalues_without_error.count(eigenvalues_without_error[ind_rho_m])
    if mult_alg_rho_m > 1:
        print(f"The spectral radius, {rho_m}, has algebraic multiplicity {mult_alg_rho_m}. "
              f"Then, is not a simple eigenvalue.")
        print(f"Therefore the matrix {matrix_name} is NOT K-irreducible.")
        return
    else:
        print(f"The spectral radius, {rho_m}, is a simple eigenvalue.")

    # We look for the eigenvalues with the same norm as the spectral radius and
    # then check if they are simple eigenvalues.

    index_norm_rho = []
    for y in range(len(norm_eigenvalues)):
        if norm_eigenvalues[y] == rho_m:
            index_norm_rho.append(y)
        else:
            continue

    for element in index_norm_rho:
        if eigenvalues_without_error.count(eigenvalues_without_error[element]) != 1:
            print(f"There is another eigenvalue, {eigenvalues_without_error[element]}, with norm {rho_m} "
                  f"and algebraic multiplicity {eigenvalues_without_error.count(eigenvalues[element])}.")
            print(f"Then, is not a simple eigenvalue. Therefore, the matrix {matrix_name} is NOT K-irreducible.")
            return
        else:
            continue
    print(f"All eigenvalues with norm {rho_m} are simples.")

    # Now check if the only eigenvector v (or -v) associated with the spectral radius is in the interior of the
    # light cone.

    U = eigenvectors.transpose()[ind_rho_m].real
    print(f"The eigenvector associated to the spectral radius is:\n{(np.array(U)[np.newaxis]).T}")
    if U[0] < 0:
        W = (-1) * U
        if ((W[3]) ** 2) + ((W[2]) ** 2) + ((W[1]) ** 2) < ((W[0]) ** 2):
            print(f"The eigenvector corresponding to the spectral radius is in the interior of the light cone.")
        else:
            print(f"The eigenvector corresponding to the spectral radius is NOT in the interior of the light cone, "
                  f"therefore {matrix_name} is NOT K-irreducible. ")
            return
    else:
        if ((U[3]) ** 2) + ((U[2]) ** 2) + ((U[1]) ** 2) < ((U[0]) ** 2):
            print(f"The eigenvector corresponding to the spectral radius is in the interior of the light light cone.")
        else:
            print(f"The eigenvector corresponding to the spectral radius is NOT in the interior of the light cone, "
                  f"therefore {matrix_name} is NOT K-irreducible. ")
            return

    # Finally, we calculate if the eigenvectors corresponding eigenvalues with norm different to the spectral radius
    # are in the light cone. If there is one, then the matrix is not K-irreducible.
    for s in range(len(eigenvalues_without_error)):
        if simplify(eigenvalues_without_error[s]).is_real == True and eigenvalues_without_error[s] != rho_m:
            U = eigenvectors.transpose()[s].real
            if U[0] < 0:
                W = (-1) * U
                if ((W[3]) ** 2) + ((W[2]) ** 2) + ((W[1]) ** 2) <= ((W[0]) ** 2):
                    print(f"The eigenvector\n{np.array(W)[np.newaxis].T} \n"
                          f"does not correspond to the spectral radius and is in the light cone. "
                          f"Therefore {matrix_name} is NOT K-irreducible.")
                    return
                else:
                    continue
            else:
                W = U
                if ((W[3]) ** 2) + ((W[2]) ** 2) + ((W[1]) ** 2) <= ((W[0]) ** 2):
                    print(f"There is an eigenvector,\n{np.array(W)[np.newaxis].T} \n"
                          f"not corresponding the spectral radius that is in the light cone. "
                          f"Therefore {matrix_name} is NOT K-irreducible.")
                    return
                else:
                    continue

        else:
            continue
    print(f"All eigenvectors not corresponding to the spectral radius are not in the light cone.\n"
          f"Therefore {matrix_name} is K-irreducible. ")


m_1 = request_matrix("M")
print(f"The input matrix M is:\n{np.array(m_1).astype(np.float64)}")
k_irreducible("M", m_1)
