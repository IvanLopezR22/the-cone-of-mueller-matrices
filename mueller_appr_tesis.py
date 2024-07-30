import numpy as np
from sympy import *
from matrix_norm import matrix_norm
from utils.request_matrix import request_matrix
from know_if_is_mueller import know_if_is_mueller

np.set_printoptions(precision=10, suppress=True, linewidth=2000)


def appr_mueller_matrix(matrix_name: str, main_matrix: Matrix) -> tuple[Matrix, str]:
    """
    Given a sympy 4x4 matrix M, this function calculates an approximation by a Mueller matrix, and we denote it by M(mu).
    The matrix M(mu)=2||M||_2*E_11+M. where ||M||_2 is the 2-norm of M and E_11 is the 4x4 matrix with the entry 11
    equal 1 and 0 elsewhere.
    The matrix M(mu) is a Mueller matrix.

    Parameters:
    -matrix_name (str): The name that the user wants to give to the matrix.
    -main_matrix (sympy 4x4 matrix): The matrix that the user wants to approximate by a Mueller matrix.

    Returns:
    Matrix: Sympy 4x4 matrix that is Mueller.
    tuple[Matrix, Str]:
        - Matrix: Sympy 4x4 matrix that is Mueller.
        - Str: M if the matrix is Mueller and M(mu) if the matrix is not Mueller.
    """
    # The matrix e_1 in the code is called E_{11} in the comments.
    e_1 = Matrix([[1, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]])
    norm_m = matrix_norm(main_matrix)
    qm_0, mi_qm_0, minimum_qm_0, proy1_m_0, mi_proy1_m_0, minimum_proy1_m_0, fig_0 = know_if_is_mueller(main_matrix)

    # In the first case we have that M is already a Mueller matrix. Then we do nothing to M.
    if minimum_qm_0 >= 0 and minimum_proy1_m_0 >= 0:
        m_appr = main_matrix
        print(f"The minimum of the functions q{matrix_name} and proy1({matrix_name}) of the matrix {matrix_name} are "
              f"{minimum_qm_0} and {minimum_proy1_m_0} respectively.")
        print(f"Therefore, the matrix {matrix_name} is a Mueller matrix. Then {matrix_name} need not to be "
              f"approximated by a Mueller matrix.")
        if __name__ == '__main__':
            print(f"The matrix {matrix_name}={matrix_name}(mu) is: \n{np.array(main_matrix).astype(np.float64)}.")
        name = str(matrix_name)

        return m_appr, name

    else:
        print(f"The minimum of the functions q{matrix_name} and proy1({matrix_name}) of the matrix {matrix_name} are "
              f"{minimum_qm_0} and {minimum_proy1_m_0} respectively.")
        print(f"Therefore, the matrix {matrix_name} is not a Mueller matrix. Then we approximate {matrix_name} to {matrix_name}(mu):")
        rang = 2 * (matrix_norm(main_matrix))  # Lo que se tiene que sumar en la primera coordenada con el primer método
        n = 1
        k = 1 / 2 * rang
        appr_new_mu = rang * e_1 + main_matrix  # Esta es la matriz con la aproximación original
        print(f"------------------------------------------------------------------------------------")
        print(np.array(appr_new_mu).astype(np.float64))
        print(f"The former matrix is mueller.")
        s = 1 / 2 * rang
        while s > 1 / (2 ** 12):
            appr_new_mu = k * e_1 + main_matrix
            print(np.array(appr_new_mu).astype(np.float64))
            print(f"The distance to the infimum is less than or equal to {s}")
            qm, mi_qm, minimum_qm, proy1_m, mi_proy1_m, minimum_proy1_m, fig = know_if_is_mueller(appr_new_mu)
            print(f"The minimum of the functions qM(mu) and proy1(M(mu)) of M(mu) are "
                  f"{minimum_qm} and {minimum_proy1_m} respectively.")
            if minimum_qm >= 0 and minimum_proy1_m >= 0:
                print(f"The former matrix is mueller")
            else:
                print(f"The former matrix is NOT Mueller")

            if minimum_qm >= 0 and minimum_proy1_m >= 0:
                x = k - (1 / (2 ** (n + 1))) * rang
            else:
                x = k + (1 / (2 ** (n + 1))) * rang
            s = (1 / (2 ** (n + 1))) * rang
            n = n + 1
            k = x
            if __name__ == '__main__':
                print(f"------------------------------------------------------------------------------------")
        qm_1, mi_qm_1, minimum_qm_1, proy1_m_1, mi_proy1_m_1, minimum_proy1_m_1, fig_1 = know_if_is_mueller(
            appr_new_mu)
        if minimum_qm_1 >= 0 and minimum_proy1_m_1 >= 0:
            m_appr = appr_new_mu
        else:
            m_appr = (1 / (2 ** 12)) * e_1 + appr_new_mu
            while not (minimum_qm_1 >= 0 and minimum_proy1_m_1 >= 0):
                m_appr = (1 / (2 ** 12)) * e_1 + m_appr
                qm_1, mi_qm_1, minimum_qm_1, proy1_m_1, mi_proy1_m_1, minimum_proy1_m_1, fig_1 = know_if_is_mueller(
                    m_appr)

        qm_1, mi_qm_1, minimum_qm_1, proy1_m_1, mi_proy1_m_1, minimum_proy1_m_1, fig_1 \
            = know_if_is_mueller(m_appr)
        print("------------------------------------------------------------------------------------")
        print(f"The approximation M(mu) of M by a Mueller matrix is:"
              f"\n{np.array(m_appr).astype(np.float64)}")
        print(f"The minimum of the functions qM(mu) and proy1(M(mu)) of M(mu) are "
              f"{minimum_qm_1} and {minimum_proy1_m_1} respectively.")
        print(f"Therefore, M(mu) is a Mueller matrix.")

        name = str(matrix_name) + "(mu)"

        return m_appr, name


if __name__ == '__main__':
    m_0 = request_matrix("M")
    print("------------------------------------------------------------------------------------")
    print(f"The input matrix M is:\n{np.array(m_0).astype(np.float64)}")
    print("------------------------------------------------------------------------------------")
    appr_mueller_matrix("M", m_0)
