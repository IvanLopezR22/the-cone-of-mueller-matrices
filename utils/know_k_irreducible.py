from sympy import *
import numpy as np
from numpy.linalg import eig
from know_if_is_mueller import know_if_is_mueller

np.set_printoptions(precision=10, suppress=True, linewidth=2000)


def k_irreducible_k_primitive(m):
    # First we need to know if the matrix is Mueller.
    qm, mi_qm, minimum_qm, proy1_m, mi_proy1_m, minimum_proy1_m, fig = know_if_is_mueller(m)
    if minimum_qm >= 0 and minimum_proy1_m >= 0:
        pass
    else:
        return False

    m_0 = np.array(m).astype(np.float64)
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

    # norm_eigenvalues is the list of the norms of the eigenvalues.
    norm_eigenvalues = (np.absolute(eigenvalues_without_error)).tolist()
    # The spectral radius ro_m is the maximum in the list norm_eigenvalue.
    rho_m = max(norm_eigenvalues)

    # Now we ask if ro_m is in the list of eigenvalues_without_error.
    if rho_m in eigenvalues_without_error:
        ind_rho_m = eigenvalues_without_error.index(rho_m)
    else:
        return False

    # We need to know if the spectral radius is a simple eigenvalue. That is, we need to know if the algebraic
    # multiplicity of the spectral radius is 1.
    mult_alg_rho_m = eigenvalues_without_error.count(eigenvalues_without_error[ind_rho_m])
    if mult_alg_rho_m > 1:
        return False
    else:
        pass
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
            return False
        else:
            continue

    # Now check if the only eigenvector v (or -v) associated with the spectral radius is in the interior of the
    # light cone.

    U = eigenvectors.transpose()[ind_rho_m].real
    if U[0] < 0:
        W = (-1) * U
        if ((W[3]) ** 2) + ((W[2]) ** 2) + ((W[1]) ** 2) < ((W[0]) ** 2):
            pass
        else:
            return False
    else:
        if ((U[3]) ** 2) + ((U[2]) ** 2) + ((U[1]) ** 2) < ((U[0]) ** 2):
            pass
        else:
            return False

    # Finally, we calculate if the eigenvectors corresponding eigenvalues with norm different to the spectral radius
    # are in the light cone. If there is one, then the matrix is not K-irreducible.
    for s in range(len(eigenvalues_without_error)):
        if simplify(eigenvalues_without_error[s]).is_real == True and eigenvalues_without_error[s] != rho_m:
            U = eigenvectors.transpose()[s].real
            if U[0] < 0:
                W = (-1) * U
                if ((W[3]) ** 2) + ((W[2]) ** 2) + ((W[1]) ** 2) <= ((W[0]) ** 2):
                    return False
                else:
                    continue
            else:
                W = U
                if ((W[3]) ** 2) + ((W[2]) ** 2) + ((W[1]) ** 2) <= ((W[0]) ** 2):
                    return False
                else:
                    continue

        else:
            continue

    return True

    #eigenvalues_norm_rho = []
    #for i in index_norm_rho:
     #   eigenvalues_norm_rho.append(eigenvalues_without_error[i])
    #if len(index_norm_rho) > 1:
     #   return True, False
    #else:
     #   return True, True
