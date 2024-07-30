import sympy as sym
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm
from sympy import *
import warnings
matplotlib.use('TkAgg')
backend = 'TkAgg'
matplotlib.use(backend)  # not detected
matplotlib.use('TkAgg', force=False)  # not detected
warnings.filterwarnings("ignore", category=RuntimeWarning)


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
    is_stokes = (first_column_m[0] ** 2) - (first_column_m[1] ** 2) - (first_column_m[2] ** 2) - (
                first_column_m[3] ** 2)

    if is_stokes >= 0 and first_column_m[0] >= 0:
        first_column_is_stokes = True
    else:
        first_column_is_stokes = False

    return qm, mi_qm, minimum_qm, first_column_is_stokes, fig
