from sympy import *
import numpy as np
from numpy.linalg import eig
from appr_invertible_matrix import make_invertible
from appr_invertible_mueller_matrix import appr_invert_mueller_matrix
from request_matrix import request_matrix
np.set_printoptions(precision=10, suppress=True, linewidth=2000)


def ecm(m):
    aw, eps_0 = make_invertible("aw", request_matrix("aw"))
    print(f"Then, the matrix aw is: \n{np.array(aw).astype(np.float64)}")
    print("------------------------------------------------------------------------------------")

    amw = request_matrix("amw")
    print('The matrix amw is: \n', np.array(amw).astype(np.float64))

    # We generate the canonical basis E_i, 1<=i<=4 of M_{4}(R) (4x4 matrices on real numbers)
    canonical_base = []
    for i in range(4):
        for j in range(4):
            def f(s, t):
                if s == i and t == j:
                    return 1
                else:
                    return 0

            k = Matrix(4, 4, f)
            canonical_base.append(k)

    aw_inv = aw.inv()

    # To represent the function H in matrix form we compute the image of the canonical basis under H.
    image_base = []
    for element in canonical_base:
        w = (m * element) - (element * aw_inv * amw)
        image_base.append(w)

    # The image of E_i under H, is the column i of H.
    h_matrix_columns = []
    for x in range(4):
        for y in range(4):
            row = []
            for element in image_base:
                row.append(element[x, y])
            h_matrix_columns.append(row)

    # matrix_h is the matrix form of H in canonical base.
    matrix_h = Matrix(h_matrix_columns)

    # matrix_m_np is the matrix H in numpy format.
    matrix_h_np = np.array(matrix_h).astype(np.float64)
    print(f"The matrix form of the function H in canonical basis is: \n{np.around(matrix_h_np, decimals=5)}")
    print("------------------------------------------------------------------------------------")

    # null_space_h is the nullspace of H.
    null_space_h = matrix_h.nullspace()

    # First consider the case when the null space of H is not trivial.
    # In this case we search for an invertible matrix in the null space of H.
    # Otherwise, we take any matrix in the nullspace of H.
    # Let W such a matrix, then we have two cases:
    # Case 1: If the matrix W is invertible, then we do nothing to W.
    # Case 2: If the matrix W is not invertible, then we force W to be invertible using make_invertible(W).
    # Note: The function make_invertible is in the file appr_mueller_matrix.
    if null_space_h:
        print('The null space of H is not trivial')
        j = 0
        while j <= (len(null_space_h) - 1) and Abs(
                round(Matrix(np.array(null_space_h[j].transpose()).astype(np.float64).reshape((4, 4))).det(),
                      15)) == 0:
            j += 1
        else:
            if j == len(null_space_h):
                a = Matrix(np.array(null_space_h[0].transpose()).astype(
                    np.float64).reshape((4, 4)))
                print(f"No invertible matrices found in the nullspace of H, so we select the matrix: W=\n"
                      f"{np.array(a).astype(np.float64)}")
                print("------------------------------------------------------------------------------------")
            else:
                a = Matrix(np.array(null_space_h[j].transpose()).astype(
                    np.float64).reshape((4, 4)))
                print(f"An invertible matrix was found in the null space of H: W=\n{np.array(a).astype(np.float64)}")
                print("------------------------------------------------------------------------------------")

        w, eps= make_invertible("W", Matrix(a))
        print(f"Then, the matrix W is: \n{np.array(w).astype(np.float64)}")
        print("------------------------------------------------------------------------------------")

    else:
        # Now we have that the null space of H is non-trivial. Then we calculate the eigenvalues of H.
        # Then there are two cases:

        # Case 1. There are no real eigenvalues of H. In this case we use the matrix H^{T}*H (the product of H transpose
        # and H). Then :
        # 1.1. We calculate the eigenvalues of H^{T}*H (all are positive values).
        # 1.2. We take the minimum eigenvalue of H^{T}*H
        # 1.3. We find an eigenvector W associated to the minimum eigenvalue of H^{T}*H.
        # 1.4. We make W invertible using the function make_invertible(W).

        # Case 2. There are real eigenvalues of H. Then:
        # 1.1. We calculate the norms of the real eigenvalues.
        # 1.2. We take the minimum of the norms of the real eigenvalues.
        # 1.3. We find an eigenvector W associated to the eigenvalue at which the minimum norm is reached.
        # 1.4. We make W invertible using the function make_invertible(W).

        print('The null space of H is trivial')
        # First we compute the eigenvalues and eigenvectors of H.

        # The 1d numpy array eigenvalues_h and the 2d numpy array eigenvectors_h
        # are the eigenvalues and eigenvectors respectively of H computed numerically.
        eigenvalues_h, eigenvectors_h = eig(matrix_h_np)

        # Numerical approximations of eigenvalues generate small imaginary parts that we need to filter.
        # For this, we use the fact that if p is a polynomial with real coefficients and a+bi is a complex root of p,
        # then a-bi is also a root of p.
        eigenvalues_h_without_error = []
        for i in eigenvalues_h:
            conj = i.conjugate()
            z = conj in eigenvalues_h
            if simplify(i).is_real == True:
                eigenvalues_h_without_error.append(i)
            elif simplify(i).is_real == False and z == False:
                eigenvalues_h_without_error.append(re(i))
            else:
                eigenvalues_h_without_error.append(i)
        print(f"The eigenvalues of H are: {eigenvalues_h_without_error}.")
        print("------------------------------------------------------------------------------------")

        # The 2d numpy array eigenvectors_h, has the eigenvectors of H arranged in the columns.
        # So for convenience we calculate the transposed matrix of eigenvectors_h
        # to have the eigenvalues of H ordered in the rows.
        eigenvectors_h_transpose = eigenvectors_h.transpose()

        # We need to know if there are real eigenvalues, so we separate the eigenvalues of H into two lists,
        # one of real eigenvalues and the other of complex eigenvalues.
        real_eigenvalues = []
        complex_eigenvalues = []
        for i in eigenvalues_h_without_error:
            if simplify(i).is_real == True:
                real_eigenvalues.append(i)
            else:
                complex_eigenvalues.append(i)
        print(f"The real eigenvalues of H are: {real_eigenvalues}.")
        print(f"The complex eigenvalues of H are: {complex_eigenvalues}.")
        print("------------------------------------------------------------------------------------")
        # Case 1:
        # There are no real eigenvalues of H.
        if not real_eigenvalues:
            # H^{T}*H is the product of the transpose of H and H.
            print('There are no real eigenvalues, therefore the approximation will be done in the matrix H^T*H.')
            ht_h = np.matmul(matrix_h_np.transpose(), matrix_h_np)

            # eigenvalues_ht_h, eigenvectors_ht_h are the eigenvalues and eigenvectors of HT^{T}H computed numerically.
            eigenvalues_ht_h, eigenvectors_ht_h = eig(ht_h)

            # The eigenvectors_ht_h_transpose matrix has the eigenvectors of H arranged in the columns.
            # Then, for convenience we calculate the transposed matrix of eigenvectors_ht_h_transpose
            # to have the eigenvalues of H ordered in the rows.
            eigenvectors_ht_h_transpose = eigenvectors_ht_h.transpose()

            print(f"The matrix H^{{T}}*H is: \n{ht_h}")
            print("------------------------------------------------------------------------------------")

            # All eigenvalues of ht_h are non-negative, but they come with a tiny imaginary part given by the numerical
            # approximation, so we take only the real part of each approximate eigenvalue.
            eigenvalues_ht_h_no_error = []
            for i in eigenvalues_ht_h:
                eigenvalues_ht_h_no_error.append(re(i))

            print(f"The eigenvalues of H^{{T}}*H are: {eigenvalues_ht_h_no_error}")
            print(f"The smallest eigenvalue of H^{{T}}*H is: {min(eigenvalues_ht_h_no_error)}")
            print("------------------------------------------------------------------------------------")

            # In this part, we find the eigenvector associated to the minimum of the norms of the eigenvalues.
            a = np.array(
                eigenvectors_ht_h_transpose[eigenvalues_ht_h_no_error.index(min(eigenvalues_ht_h_no_error))]).astype(
                np.float64).reshape((4, 4))
            print(f"An eigenvector W associated to the eigenvalue {min(eigenvalues_ht_h_no_error)} is: \n{a}")
            print("------------------------------------------------------------------------------------")

            w, eps = make_invertible("W", Matrix(a))
            print(f"Then, the matrix W is: \n{np.array(w).astype(np.float64)}")
            print("------------------------------------------------------------------------------------")

        else:
            # Case 2:
            # There are real eigenvalues of H.
            # We need to find the real eigenvalue with the smallest norm.
            abs_real_eigenvalues = []
            for i in range(len(real_eigenvalues)):
                abs_real_eigenvalues.append(Abs(real_eigenvalues[i]))

            print(f"The norms of the real eigenvalues of H are: {abs_real_eigenvalues}")
            print(f"The minimum of the norms of the eigenvalues is: {min(abs_real_eigenvalues)}")
            print("------------------------------------------------------------------------------------")

            ind_min_real_eigenvalues = eigenvalues_h_without_error.index(
                real_eigenvalues[abs_real_eigenvalues.index(min(abs_real_eigenvalues))])
            print(f"An eigenvalue associated with the minimum norm is:"
                  f"\n{eigenvalues_h_without_error[ind_min_real_eigenvalues]}")

            # The matrix called "a" in the code is the matrix A associated to the minimum norm.
            a = np.array(eigenvectors_h_transpose[ind_min_real_eigenvalues]).reshape((4, 4))
            print(f"An eigenvector W associated to the eigenvalue "
                  f"{eigenvalues_h_without_error[ind_min_real_eigenvalues]} is:\n{a}")
            print("------------------------------------------------------------------------------------")
            w, eps = make_invertible("W", Matrix(a))
            print(f"Then, the matrix W is: \n{np.array(w).astype(np.float64)}")
            print("------------------------------------------------------------------------------------")

    # The last step is to calculate a new M:
    # 1.1. Let new_M=W*(aw^{-1})*amw*(W^{-1})
    # 1.2. We approximate new_M by an invertible Mueller matrix using the function
    # new_M(mu-inv)=appr_invert_mueller_matrix(M_{0})
    # Then, new_M(mu-inv) is our final matrix, an invertible Mueller matrix.

    new_m = w*aw_inv*amw*(w.inv())
    print()
    print(f"The new matrix M, denoted by New_M, is the matrix defined by new_M=W*(aw^{-1})*amw*(W^{-1})="
          f"\n{np.array(new_m).astype(np.float64)}")
    print("------------------------------------------------------------------------------------")
    print(f"Following the notation of the comments in the appr_invertible_mueller_matrix and "
          f"appr_by_mueller_matrix files, we call New_M(mu-inv) to the approximation ")
    print(f"of New_M by an invertible Mueller matrix. This is the final matrix obtained by the ECM.")
    print("------------------------------------------------------------------------------------")

    m_appr = appr_invert_mueller_matrix("New_M", Matrix(np.around(np.array(new_m).astype(np.float64), decimals=8)))

    return m_appr


new_m_0 = request_matrix("M")
print('The matrix M is: \n', np.array(new_m_0).astype(np.float64))
ecm(new_m_0)
