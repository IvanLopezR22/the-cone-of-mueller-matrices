from sympy import *
from sympy import Abs
import numpy as np
from numpy.linalg import eig
import sympy as sym
import matplotlib.pyplot as plt
import warnings

warnings.filterwarnings("ignore", category=RuntimeWarning)

np.set_printoptions(precision=15, suppress=True, linewidth=2000)


def request_matrix(matrix_name):
    print(f"Enter the values of the matrix {matrix_name}:")
    matrix_m = []
    for raw_position in range(4):
        raw = []
        for element in range(4):
            while True:
                try:
                    value = Rational(str(float((input(f"Enter the element {raw_position + 1}, {element + 1}: ")))))
                    if raw_position == 0 and element == 0 and value == 0:
                        print("The value at coordinate (1, 1) must be non-zero. Please try again.")
                        continue

                    raw.append(value)
                except ValueError:
                    print("A non-numeric character was entered, please input a decimal real number.")
                    continue

                else:
                    break

        matrix_m.append(raw)

    return Matrix(matrix_m)


def matrix_norm(main_matrix):
    mult_transpose_main = main_matrix.T * main_matrix
    eigenvalues_h, eigenvectors_h = eig(np.array(mult_transpose_main).astype(np.float64))
    norm_eigenvalues_mtm = []
    for i in range(len(eigenvalues_h)):
        norm_eigenvalues_mtm.append(Abs(N(re(eigenvalues_h[i]), 10)))
    s = sqrt(max(norm_eigenvalues_mtm))
    return s


def van_der_mee_theorem(main_matrix):
    G = Matrix([[1, 0, 0, 0], [0, -1, 0, 0], [0, 0, -1, 0], [0, 0, 0, -1]])
    x, y = sym.symbols('x,y')
    s = main_matrix * Matrix([1, x, y, sqrt(1 - x ** 2 - y ** 2)])
    t = G * (main_matrix * Matrix([1, x, y, sqrt(1 - x ** 2 - y ** 2)]))
    qm = s.dot(t)
    proy1_m = (main_matrix * Matrix([1, x, y, sqrt(1 - x ** 2 - y ** 2)]))[0]
    z_lamb_1 = sym.lambdify((x, y), qm)
    z_lamb_2 = sym.lambdify((x, y), proy1_m)
    a = 1
    b = 1
    xdata = np.linspace(-(b + .25), (b + .25), 1001)
    ydata = np.linspace(-(b + .25), (b + .25), 1001)
    X, Y = np.meshgrid(xdata, ydata)
    fig = plt.figure(figsize=(3.5 * 3.13, 1.5 * 3.13))
    R = (X ** 2 + Y ** 2 > a ** 2)
    expr_1 = sympify(qm)
    symbols_qm = expr_1.free_symbols
    expr_2 = sympify(proy1_m)
    symbols_m_1 = expr_2.free_symbols
    if not symbols_qm:
        K_1 = np.around(np.full((1001, 1001), qm).astype(np.float64), decimals=12)
    else:
        K_1 = np.around(z_lamb_1(X, Y), decimals=12)

    if not symbols_m_1:
        K_2 = np.around(np.full((1001, 1001), proy1_m).astype(np.float64), decimals=12)
    else:
        K_2 = np.around(z_lamb_2(X, Y), decimals=12)
    z_masked_1 = np.ma.masked_where(R, K_1)
    z_masked_2 = np.ma.masked_where(R, K_2)
    xmin_qm, ymin_qm = np.unravel_index(np.argmin(z_masked_1), z_masked_1.shape)
    mi_qm = (X[xmin_qm, ymin_qm], Y[xmin_qm, ymin_qm], z_masked_1.min())
    minimum_qm = z_masked_1.min()
    xmin_proy1_m, ymin_proy1_m = np.unravel_index(np.argmin(z_masked_2), z_masked_2.shape)
    mi_proy1_m = (X[xmin_proy1_m, ymin_proy1_m], Y[xmin_proy1_m, ymin_proy1_m], z_masked_2.min())
    minimum_proy1_m = z_masked_2.min()
    return qm, mi_qm, minimum_qm, proy1_m, mi_proy1_m, minimum_proy1_m, fig


def appr_mueller_matrix(w, m):
    e_1 = Matrix([[1, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]])
    norm_m = matrix_norm(m)
    qm, mi_qm, minimum_qm, proy1_m, mi_proy1_m, minimum_proy1_m, fig = van_der_mee_theorem(m)
    if minimum_qm >= 0 and minimum_proy1_m >= 0:
        m_appr = m
        print(f"The minimum of the functions q{w} and proy1({w}) of the matrix {w} are "
              f"{minimum_qm} and {minimum_proy1_m} respectively.")
        print(f"Therefore, the matrix {w} is a Mueller matrix. Then {w} need not to be "
              f"approximated by a Mueller matrix.")
        name = str(w)
        return m_appr, name

    else:
        m_appr = ((2 * norm_m) * e_1) + m
        qm_1, mi_qm_1, minimum_qm_1, proy1_m_1, mi_proy1_m_1, minimum_proy1_m_1, fig_1 \
            = van_der_mee_theorem(m_appr)
        print(f"The minimum of the functions q{w} and proy1({w}) of the matrix {w} are "
              f"{minimum_qm} and {minimum_proy1_m} respectively.")
        print(f"Therefore, the matrix {w} is not a Mueller matrix. Then we approximate {w} to {w}(mu).")

        name = str(w) + "(mu)"

        return m_appr, name


def make_invertible(name, main_matrix):
    ident = Matrix([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    if round(main_matrix.det(), 12) != 0:
        invertible_main_matrix = main_matrix
        print(f"The determinant of the matrix {name} is {round(main_matrix.det(), 12)}. Therefore {name} "
              f"is an invertible matrix. ")
    else:
        eigen_m = main_matrix.eigenvects()
        eigenvalues = []
        for element in eigen_m:
            eigenvalues.append((element[0]))
        norm_eigenvalues_m = []
        for i in range(len(eigen_m)):
            norm_eigenvalues_m.append(round(Abs(N(eigen_m[i][0], 12)), 12))
        norms_no_zero = [x for x in norm_eigenvalues_m if x != 0]
        if not norms_no_zero:
            norms_no_zero.append(2)
        eps = min([1 / 100, (1 / 2) * min(norms_no_zero)])
        invertible_main_matrix = (eps * ident) + main_matrix
        print(f"The determinant of the matrix {name} is {round(main_matrix.det(), 12)}. "
              f"Therefore {name} is not an invertible matrix. ")
        print("------------------------------------------------------------------------------------")

    return invertible_main_matrix


def appr_invert_mueller_matrix(w, main_matrix):
    m_mue_appr, name = appr_mueller_matrix(w, main_matrix)
    m_mue_inv = make_invertible(name, m_mue_appr)
    qm, mi_qm, minimum_qm, proy1_m, mi_proy1_m, minimum_proy1_m, fig = van_der_mee_theorem(m_mue_inv)
    print(f"Due to the previous data, the approximation {w}(mu-inv) of {w} by an invertible Mueller matrix is: "
          f"\n{np.array(m_mue_inv).astype(np.float64)}")
    print(f"The minimum of the functions q{w}(mu-inv) and proy1({w}(mu-inv)) of {w}(mu-inv) are "
          f"{minimum_qm} and {minimum_proy1_m} respectively.")
    print(f"The determinant of {w}(mu-inv) is {m_mue_inv.det()}. Therefore, {w}(mu-inv) is "
          f"an invertible Mueller matrix.")
    return m_mue_inv


def ecm(m):
    aw = make_invertible("aw", request_matrix("aw"))
    print(f"Then, the matrix aw is: \n{np.array(aw).astype(np.float64)}")
    print("------------------------------------------------------------------------------------")
    amw = request_matrix("amw")
    print('The matrix amw is: \n', np.array(amw).astype(np.float64))
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
    image_base = []
    for element in canonical_base:
        w = (m * element) - (element * aw_inv * amw)
        image_base.append(w)
    h_matrix_columns = []
    for x in range(4):
        for y in range(4):
            row = []
            for element in image_base:
                row.append(element[x, y])
            h_matrix_columns.append(row)
    matrix_h = Matrix(h_matrix_columns)
    matrix_h_np = np.array(matrix_h).astype(np.float64)
    print(f"The matrix form of the function H in canonical basis is: \n{np.around(matrix_h_np, decimals=5)}")
    print("------------------------------------------------------------------------------------")
    null_space_h = matrix_h.nullspace()
    if null_space_h:
        print('The null space of H is non-trivial.')
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
        w = make_invertible("W", Matrix(a))
        print(f"Then, the matrix W is: \n{np.array(w).astype(np.float64)}")
        print("------------------------------------------------------------------------------------")

    else:
        print('The null space of H is trivial.')
        eigenvalues_h, eigenvectors_h = eig(matrix_h_np)
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
        eigenvectors_h_transpose = eigenvectors_h.transpose()
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
        if not real_eigenvalues:
            print('There are no real eigenvalues, therefore the approximation will be done in the matrix H^T*H.')
            ht_h = np.matmul(matrix_h_np.transpose(), matrix_h_np)
            eigenvalues_ht_h, eigenvectors_ht_h = eig(ht_h)
            eigenvectors_ht_h_transpose = eigenvectors_ht_h.transpose()
            print(f"The matrix H^{{T}}*H is: \n{ht_h}")
            print("------------------------------------------------------------------------------------")
            eigenvalues_ht_h_no_error = []
            for i in eigenvalues_ht_h:
                eigenvalues_ht_h_no_error.append(re(i))
            print(f"The eigenvalues of H^{{T}}*H are: {eigenvalues_ht_h_no_error}")
            print(f"The smallest eigenvalue of H^{{T}}*H is: {min(eigenvalues_ht_h_no_error)}")
            print("------------------------------------------------------------------------------------")
            a = np.array(
                eigenvectors_ht_h_transpose[eigenvalues_ht_h_no_error.index(min(eigenvalues_ht_h_no_error))]).astype(
                np.float64).reshape((4, 4))
            print(f"An eigenvector W associated to the eigenvalue {min(eigenvalues_ht_h_no_error)} is: \n{a}")
            print("------------------------------------------------------------------------------------")
            w = make_invertible("W", Matrix(a))
            print(f"Then, the matrix W is: \n{np.array(w).astype(np.float64)}")
            print("------------------------------------------------------------------------------------")
        else:
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
            a = np.array(eigenvectors_h_transpose[ind_min_real_eigenvalues]).reshape((4, 4))
            print(f"An eigenvector W associated to the eigenvalue "
                  f"{eigenvalues_h_without_error[ind_min_real_eigenvalues]} is:\n{a}")
            print("------------------------------------------------------------------------------------")
            w = make_invertible("W", Matrix(a))
            print(f"Then, the matrix W is: \n{np.array(w).astype(np.float64)}")
            print("------------------------------------------------------------------------------------")
    new_m = w * aw_inv * amw * (w.inv())
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