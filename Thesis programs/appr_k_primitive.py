import numpy as np
from sympy import *
from matrix_norm import matrix_norm
from know_if_is_mueller import know_if_is_mueller
from appr_mueller_matrix import appr_mueller_matrix
from request_matrix import request_matrix


def appr_k_primitive(matrix_name: str, main_matrix: Matrix) -> Matrix:
    e_1 = Matrix([[1, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]])
    # In this function, we make an approximation of M by a primitive matrix, denoted as M(prim).
    # Case 1. The introduced matrix M is K-primitive, then M(prim)=M
    # Case 2. The introduced matrix M is a Mueller matrix but not K-primitive, then M(prim)=(1/100)E_{11}+M.
    # Case 3. The introduced matrix M is not a Mueller matrix, then M(prim)=(1/100)E_{11}+M(mu).
    norm_m = matrix_norm(main_matrix)
    # The matrix is K-primitive if and only if the minimums of the functions qM and proj1(M) are positive.
    qm, mi_qm, minimum_qm, first_column_is_stokes_m, fig = know_if_is_mueller(main_matrix)

    # This is the case in which the matrix M is K-primitive, then we leave M unchanged.
    if minimum_qm > 0 and first_column_is_stokes_m == True:
        m_primitive = main_matrix
        print(f"The matrix {matrix_name} is K-primitive. Then the approximation is the same {main_matrix}:"
              f" \n{np.array(main_matrix).astype(np.float64)}")

        return m_primitive

    # The second case is when the matrix M is a Mueller matrix and is not K-primitive.
    # Then we have an approximation M(prim) that we call the approximation of M by a K-primitive matrix.
    elif minimum_qm > 0 and first_column_is_stokes_m == True:
        m_primitive = (1 / 100 * e_1) + main_matrix

        print(f"The matrix {matrix_name} is a Mueller matrix but is not K-primitive.")
        print(f"The approximation {matrix_name}(prim) of {matrix_name} by a K-primitive Mueller matrix is: "
              f"\n{np.array(m_primitive).astype(np.float64)}")
        qm, mi_qm, minimum_qm, first_column_is_stokes_m, fig = know_if_is_mueller(m_primitive)
        print(f"The minimum of the function q{matrix_name}(prim) of {matrix_name}(prim) is {minimum_qm}.")
        print(f"The first column of the matrix {matrix_name} is a Stokes vector: {first_column_is_stokes_m}.")
        print(f"Therefore {matrix_name}(prim) is K-primitive.")
        print(f"-----------------------------------------------------------------------------")

        return m_primitive

    else:
        # In the last case the matrix M is not a Mueller matrix, therefore is not K-primitive.
        # Then first we approximate M by a Mueller matrix and then we appy the K-primitive
        # approximation to M(mu).
        print(f"The minimum of the function q{matrix_name} of the matrix {matrix_name} is {minimum_qm}")
        print(f"The first column of the matrix {matrix_name} is a Stokes vector: {first_column_is_stokes_m}.")
        print(f"Therefore, M is not a Mueller matrix and in consequence is not K-primitive.")
        print(f"----------------------------------------------------------------------------")
        print(f"First we approximate {matrix_name} by a Mueller matrix:")
        m_mueller, name = appr_mueller_matrix(matrix_name, main_matrix)
        print(f"----------------------------------------------------------------------------")
        print(f"Now we approximate to a K-primitive matrix.")
        m_primitive = (m_mueller) + (1 / 100 * e_1)
        print(f"The approximation {matrix_name}(prim) of {matrix_name} by a K-primitive Mueller matrix is: "
              f"\n{np.array(m_primitive).astype(np.float64)}")
        qm, mi_qm, minimum_qm, first_column_is_stokes_m, fig = know_if_is_mueller(m_primitive)
        print(f"The minimum of the function q{matrix_name}(prim) of {matrix_name}(prim) is {minimum_qm}.")
        print(f"The first column of the matrix {matrix_name} is a Stokes vector: {first_column_is_stokes_m}.")
        print(f"Therefore {matrix_name}(prim) is K-primitive.")
        print(f"-----------------------------------------------------------------------------")

        return m_primitive


if __name__ == '__main__':
    m_0 = request_matrix("M")
    print("------------------------------------------------------------------------------------")
    print(f"The input matrix M is:\n{np.array(m_0).astype(np.float64)}")
    print("------------------------------------------------------------------------------------")
    appr_k_primitive("M", m_0)
