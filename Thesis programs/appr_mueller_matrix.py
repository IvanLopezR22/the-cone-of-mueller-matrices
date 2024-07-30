import numpy as np
from sympy import *
from matrix_norm import matrix_norm
from request_matrix import request_matrix
from know_if_is_mueller import know_if_is_mueller

np.set_printoptions(precision=10, suppress=True, linewidth=2000)


def appr_mueller_matrix(matrix_name: str, main_matrix: Matrix) -> tuple[Matrix, str]:
    """
    Given a sympy 4x4 matrix M, this function calculates an approximation by a Mueller matrix, and we denote it by M(mu).
    The matrix M(mu)=2||M||_2*E_11+M. where ||M||_2 is the 2-norm of M and E_11 is the 4x4 matrix with the entry 11
    equal 1 and 0 elsewhere.
    In the interval [0,2||M||_2] there exist a minimum number b such that b*E_11+M is a Mueller matrix.
    We calculate a number c such that c*E_11+M is a Mueller matrix and |c-d |<1/2^6.
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
    qm_0, mi_qm_0, minimum_qm_0, first_column_is_stokes_0, fig_0 = know_if_is_mueller(main_matrix)

    # In the first case we have that M is already a Mueller matrix. Then we do nothing to M.
    if minimum_qm_0 >= 0 and first_column_is_stokes_0 == True:
        m_appr = main_matrix
        print(f"The minimum of the function q{matrix_name} of the matrix {matrix_name} is {minimum_qm_0}.")
        print(f"The first column of the matrix {matrix_name} is a Stokes vector: {first_column_is_stokes_0}.")
        print(f"Therefore, the matrix {matrix_name} is a Mueller matrix. Then {matrix_name} need not to be "
              f"approximated by a Mueller matrix.")
        if __name__ == '__main__':
            print(f"The matrix {matrix_name}={matrix_name}(mu) is: \n{np.array(main_matrix).astype(np.float64)}.")
        name = str(matrix_name)

        return m_appr, name

    else:
        print(f"The minimum of the function q{matrix_name} of the matrix {matrix_name} is {minimum_qm_0}.")
        print(f"The first column of the matrix {matrix_name} is a Stokes vector: {first_column_is_stokes_0}.")
        print(f"Therefore, the matrix {matrix_name} is not a Mueller matrix. "
              f"Then we approximate {matrix_name} to {matrix_name}(mu):")
        rang = 2 * (matrix_norm(main_matrix))
        n = 1
        k = 1 / 2 * rang
        appr_new_mu = rang * e_1 + main_matrix
        print(f"------------------------------------------------------------------------------------")
        print(np.array(appr_new_mu).astype(np.float64))
        print(f"The former matrix is mueller.")
        print(f"------------------------------------------------------------------------------------")
        s = 1 / 2 * rang
        while s > 1 / (2 ** 12):
            appr_new_mu = k * e_1 + main_matrix
            print(np.array(appr_new_mu).astype(np.float64))
            print(f"The distance to {matrix_name}(mu) is less than or equal to {s}")
            qm, mi_qm, minimum_qm, first_column_is_stokes, fig = know_if_is_mueller(appr_new_mu)
            print(f"The minimum of the function q({k}e_1+{matrix_name}) of ({k}e_1+{matrix_name}) is: {minimum_qm}.")
            print(f"The first column of the matrix {matrix_name} is a Stokes vector: {first_column_is_stokes}.")
            if minimum_qm >= 0 and first_column_is_stokes == True:
                print(f"The former matrix is mueller")
            else:
                print(f"The former matrix is NOT Mueller")

            if minimum_qm >= 0 and first_column_is_stokes == True:
                x = k - (1 / (2 ** (n + 1))) * rang
            else:
                x = k + (1 / (2 ** (n + 1))) * rang
            s = (1 / (2 ** (n + 1))) * rang
            n = n + 1
            k = x
            print(f"------------------------------------------------------------------------------------")
        qm_1, mi_qm_1, minimum_qm_1, first_column_is_stokes_1, fig_1 = know_if_is_mueller(
            appr_new_mu)
        if minimum_qm_1 >= 0 and first_column_is_stokes_1 == True:
            m_appr = appr_new_mu
        else:
            m_appr = (1 / (2 ** 12)) * e_1 + appr_new_mu
            while not (minimum_qm_1 >= 0 and first_column_is_stokes_1 == True):
                m_appr = (1 / (2 ** 12)) * e_1 + m_appr
                qm_1, mi_qm_1, minimum_qm_1, first_column_is_stokes_1, fig_1 = know_if_is_mueller(
                    m_appr)

        qm_1, mi_qm_1, minimum_qm_1, first_column_is_stokes_m_1, fig_1 \
            = know_if_is_mueller(m_appr)
        print("------------------------------------------------------------------------------------")
        print(f"The approximation {matrix_name}(mu) of M by a Mueller matrix is:"
              f"\n{np.array(m_appr).astype(np.float64)}")
        print(f"The minimum of the function {matrix_name}(mu) of {matrix_name}(mu) is {minimum_qm_1}")
        print(f"The first column of the matrix {matrix_name} is a Stokes vector: {first_column_is_stokes_m_1}.")
        print(f"Therefore, {matrix_name}(mu) is a Mueller matrix.")

        name = str(matrix_name) + "(mu)"

        return m_appr, name


if __name__ == '__main__':
    m_0 = request_matrix("M")
    print("------------------------------------------------------------------------------------")
    print(f"The input matrix M is:\n{np.array(m_0).astype(np.float64)}")
    print("------------------------------------------------------------------------------------")
    appr_mueller_matrix("M", m_0)
