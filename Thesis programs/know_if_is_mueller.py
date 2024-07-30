import sympy as sym
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from sympy import *
from request_matrix import request_matrix
import warnings

warnings.filterwarnings("ignore", category=RuntimeWarning)
np.set_printoptions(precision=10, suppress=True, linewidth=2000)


def know_if_is_mueller(main_matrix):
    """
    This function uses the van der Mee's quadratic form to calculate if a 4x4 matrix is a Mueller matrix, that is,
    calculate the function q_M and proy_1(M) of the van der Mee's quadratic form and its respective minimum.
    If both minimum are non-negative, the matrix M is Mueller, otherwise the matrix is not Mueller.

    Parameters:
    -main_matrix (sympy 4x4 matrix): The matrix that the user wants to know if is Mueller.

    Returns:
    Matrix: Sympy 4x4 matrix that is Mueller.
    tuple[Add, tuple, float, Add, tuple, float, Figure]:
        - The expression of the function q_M in terms of the variables x and y.
        - The point of R^{3} where the minimum of the function q_M is reached.
        - The minimum of the function q_M.
        - The expression of the function proy_1(M) in terms of the variables x and y.
        - The point of R^{3} where the minimum of the function proy_1(M) is reached.
        - The minimum of the function proy_1(M).
        - A Matplotlib Figure object representing the plot of the functions q_M and proy_1(M).

    """
    G = Matrix([[1, 0, 0, 0], [0, -1, 0, 0], [0, 0, -1, 0], [0, 0, 0, -1]])
    x, y = sym.symbols('x,y')
    s = main_matrix * Matrix([1, x, y, sqrt(1 - x ** 2 - y ** 2)])
    t = G * (main_matrix * Matrix([1, x, y, sqrt(1 - x ** 2 - y ** 2)]))

    # For a 4x4 matrix M, the function qM is the quadratic form of the Van Der Mee Theorem and proy1(M)
    # the projection on the first coordinate, which in the code are called qm and proy1_m respectively.
    # A matrix M in M_{4}(R) is a Mueller matrix if and only if qM>=0 and proy1(M)>=0.
    qm = s.dot(t)
    z_lamb_1 = sym.lambdify((x, y), qm)
    a = 1
    b = 1
    xdata = np.linspace(-(b + .25), (b + .25), 1001)
    ydata = np.linspace(-(b + .25), (b + .25), 1001)
    X, Y = np.meshgrid(xdata, ydata)
    fig = plt.figure(figsize=(3.5 * 3.13, 2 * 3.13))
    ax_1 = fig.add_subplot(111, projection='3d')
    R = (X ** 2 + Y ** 2 > a ** 2)
    expr_1 = sympify(qm)
    symbols_qm = expr_1.free_symbols
    if not symbols_qm:
        K_1 = np.around(np.full((1001, 1001), qm).astype(np.float64), decimals=12)
    else:
        K_1 = np.around(z_lamb_1(X, Y), decimals=12)

    z_masked_1 = np.ma.masked_where(R, K_1)
    surf_1 = ax_1.plot_surface(X, Y, z_masked_1, cmap=cm.winter, linewidth=0, antialiased=False)
    fig.colorbar(surf_1, shrink=0.5, aspect=15, location='left')
    ax_1.set_title(r'Graph of the function $q_{M}:\mathbb{D}\rightarrow \mathbb{R}$')
    ax_1.set_xlabel('x-axis')
    ax_1.set_ylabel('y-axis')
    ax_1.set_zlabel('z-axis')

    fig.suptitle("van der Mee's quadratic form", fontsize=16)
    xmin_qm, ymin_qm = np.unravel_index(np.argmin(z_masked_1), z_masked_1.shape)
    mi_qm = (X[xmin_qm, ymin_qm], Y[xmin_qm, ymin_qm], z_masked_1.min())
    ax_1.scatter(X[xmin_qm, ymin_qm], Y[xmin_qm, ymin_qm], z_masked_1.min(), c='red', s=10)
    minimum_qm = z_masked_1.min()

    first_column_m = main_matrix.col(0)
    is_stokes = (first_column_m[0]**2)-(first_column_m[1]**2)-(first_column_m[2]**2)-(first_column_m[3]**2)

    if is_stokes>=0 and first_column_m[0] >= 0:
        first_column_is_stokes = True
    else:
        first_column_is_stokes = False

    return qm, mi_qm, minimum_qm, first_column_is_stokes, fig


if __name__ == '__main__':
    m_0 = request_matrix("M")
    print(f"The input matrix M is:\n{np.array(m_0).astype(np.float64)}")
    qm_0, mi_qm_0, minimum_qm_0, first_column_is_stokes, fig_0 = know_if_is_mueller(m_0)
    if minimum_qm_0 >= 0 and first_column_is_stokes == True:
        result = (f"The first column of the matrix M is a Stokes vector: {first_column_is_stokes}.\n"
                  f"The minimum of the function qM of the matrix M is {minimum_qm_0}.\n"
                  f"Therefore M is a Mueller matrix.")
        fig_0.text(0.5, 0.02, result, ha='center', fontsize=12)
    else:
        result = (f"The first column of the matrix M is a Stokes vector: {first_column_is_stokes}.\n"
                  f"The minimum of the function qM of the matrix M is {minimum_qm_0}.\n"
                  f"Therefore M is NOT a Mueller matrix.")
        fig_0.text(0.5, 0.02, result, ha='center', fontsize=12)
    print(f"-------------------------------------------------------------------------------")
    print(f"The function qm of the matrix M is:\nqm= {qm_0}.")
    print(f"A minimum point on the graph of qm is: {mi_qm_0}.")
    print(f"-------------------------------------------------------------------------------")
    print(f"The minimum value of the function qm is {minimum_qm_0}.")
    print(f"The first column of M is a Stokes vector: {first_column_is_stokes}.")
    print(f"-------------------------------------------------------------------------------")

    min_stock = Matrix([1, mi_qm_0[0], mi_qm_0[1], sqrt(round(1 - mi_qm_0[0] ** 2 - mi_qm_0[1] ** 2, 10))])
    print(
        f"The minimum of the function qM in the light cone is the vector v defined by : \n{np.array(min_stock).astype(np.float64)}.")
    mat_1 = m_0 * min_stock
    print(f"The result of the multiplication of M and v is:\n {np.array(mat_1).astype(np.float64)} ")
    # print(f"{np.array(main_matrix).astype(np.float64)}{np.array(mat_0).astype(np.float64)}={np.array(mat_1).astype(np.float64)}")
    ineq = mat_1[1] ** 2 + mat_1[2] ** 2 + mat_1[3] ** 2
    if ineq > (mat_1[0] ** 2):
        print(f"{ineq}={mat_1[1]} ** {2} + {mat_1[2]} ** {2}+ {mat_1[3]} ** {2}> {mat_1[0]} ** {2} = {mat_1[0] ** 2}")
        print(f"Then v is not a Stock's vector")
    else:
        print(f"{ineq}={mat_1[1] ** 2} + {mat_1[2] ** 2}+ {mat_1[3] ** 2}<= {mat_1[0]} = {mat_1[0] ** 2}")
        print(f"Then v is a Stock's vector")
    if minimum_qm_0 < 0 or first_column_is_stokes < 0:
        print('The matrix M is NOT a Mueller matrix.')
    else:
        print('The matrix M is a Mueller matrix.')

    plt.show()
