import sympy as sym
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm
from sympy import *
matplotlib.use('TkAgg')
backend = 'TkAgg'
matplotlib.use(backend)  # not detected

matplotlib.use('TkAgg', force=False)  # not detected


def qm_graphic(main_matrix):
    G = Matrix([[1, 0, 0, 0], [0, -1, 0, 0], [0, 0, -1, 0], [0, 0, 0, -1]])
    x, y = sym.symbols('x,y')
    s = main_matrix * Matrix([1, x, y, sqrt(1 - x ** 2 - y ** 2)])
    t = G * (main_matrix * Matrix([1, x, y, sqrt(1 - x ** 2 - y ** 2)]))

    # qm is the quadratic form function of the Van Der Mee Theorem.
    qm = s.dot(t)
    print(f"The function qm of the matrix M is:\n qm= {qm}.")

    # Next we plot the function qm.
    z_lamb = sym.lambdify((x, y), qm)
    a = 1
    b = 1
    xdata = np.linspace(-(b + .25), (b + .25), 1001)
    ydata = np.linspace(-(b + .25), (b + .25), 1001)
    X, Y = np.meshgrid(xdata, ydata)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    R = (X ** 2 + Y ** 2 > a ** 2)
    expr = sympify(qm)
    symbols_qm = expr.free_symbols
    if not symbols_qm:
        K = np.full((1001, 1001), qm).astype(np.float64)
    else:
        K = z_lamb(X, Y)
    z_masked = np.ma.masked_where(R, K)
    surf = ax.plot_surface(X, Y, z_masked, cmap=cm.winter,
                           linewidth=0, antialiased=False)
    fig.colorbar(surf, shrink=0.5, aspect=5)
    return z_masked, ax, qm, X, Y
