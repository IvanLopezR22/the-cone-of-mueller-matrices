from sympy import *


def request_matrix(matrix_name: str) -> Matrix:
    """
    This function allow the user to enter a 4x4 Sympy matrix or Matrix.

    Parameters:
    -matrix_name (string): The name of the matrix.

    Returns:
    A Sympy 4x4 matrix.

    Raises:
    ValueError: If the entry (1,1) of the matrix is zero.
    ValueError: A non-numeric character was entered.
    """
    print(f"Enter the values of the matrix {matrix_name}:")
    matrix_m = []
    for raw_position in range(4):
        raw = []
        for element in range(4):
            while True:
                try:
                    value = float((input(f"Enter the element {raw_position + 1}, {element + 1}: ")))
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