import numpy as np
from sympy import *

np.set_printoptions(precision=10, suppress=True, linewidth=2000)


def make_invertible(name: str, main_matrix: Matrix):
    """ Let M a 4x4 matrix with entries in the real numbers. This function finds and epsilon>=0
        such that M(inv,epsilon):=M+epsilon*Id_4 is invertible, where Id_4 is the identity matrix of dimension 4.

        Parameters:
        -name (string): The name of the matrix.
        -main_matrix (sympy 4x4 Matrix): The matrix you want to approximate by an invertible matrix.

        Returns:
        A sympy 4x4 matrix, that is a modification of main_matrix, such that is invertible.
    """

    ident = Matrix([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])

    # If M is a non-invertible matrix, then we call M(inv)=epsilon*Id+M.

    # Case 1: The matrix M is invertible, then we do nothing to M.
    det_main_matrix = round(main_matrix.det(), 12)
    if det_main_matrix != 0:
        invertible_main_matrix = main_matrix
        print(f"The determinant of the matrix {name} is {round(main_matrix.det(), 12)}. Therefore {name} "
              f"is an invertible matrix. ")
        if __name__ == '__main__':
            print(f"Then, {name}(inv,0)={name} is: \n{np.array(invertible_main_matrix).astype(np.float64)}")
        print("------------------------------------------------------------------------------------")

    # Case 2: The matrix M is not invertible. Then we calculate the epsilon mentioned above.
    else:
        # eigen_m is the list of eigenvalues, algebraic multiplicity and eigenvectors of the matrix M.
        eigen_m = main_matrix.eigenvects()
        eigenvalues = []
        for element in eigen_m:
            eigenvalues.append((element[0]))

        # Now we calculate the norm of the eigenvalues.
        norm_eigenvalues_m = []
        for i in range(len(eigen_m)):
            norm_eigenvalues_m.append(round(Abs(N(eigen_m[i][0], 10)), 15))

        # Here we filter the zero eigenvalues.
        norms_no_zero = [x for x in norm_eigenvalues_m if x != 0]
        if not norms_no_zero:
            norms_no_zero.append(2)
        # Using the list of norms of eigenvalues we find our epsilon.
        eps = min([1 / 100, (1 / 2) * min(norms_no_zero)])
        # Finally we add epsilon*ident to our input matrix.
        invertible_main_matrix = (eps * ident) + main_matrix
        print(f"The determinant of the matrix {name} is {round(main_matrix.det(), 12)}. "
              f"Therefore {name} is not an invertible matrix. ")
        print("------------------------------------------------------------------------------------")
        if __name__ == '__main__':
            print(f"The matrix {name}(inv,{eps}) is: \n{np.array(invertible_main_matrix).astype(np.float64)}")
            print(f"The matrix {name}(inv,{eps}) is invertible, because det({name}(inv,{eps})="
                  f"{invertible_main_matrix.det()}")

    return invertible_main_matrix, det_main_matrix
