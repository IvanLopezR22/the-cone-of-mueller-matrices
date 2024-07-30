import numpy as np
from matrix_norm import *
from sympy import *
from utils.request_matrix import request_matrix_interface
from utils.format_print_matrix import format_print_matrix
from know_if_is_mueller import know_if_is_mueller
from appr_invertible_mueller_matrix import appr_invert_mueller_matrix
from ecm import ecm
from k_irreducible import k_irreducible
from k_primitive import k_primitive
from appr_invertible_matrix import make_invertible
from appr_mueller_matrix import appr_mueller_matrix
from appr_k_primitive import appr_k_primitive
import matplotlib.pyplot as plt
import ctypes


LF_FACESIZE = 32
STD_OUTPUT_HANDLE = -11


class COORD(ctypes.Structure):
    _fields_ = [("X", ctypes.c_short), ("Y", ctypes.c_short)]


class CONSOLE_FONT_INFOEX(ctypes.Structure):
    _fields_ = [("cbSize", ctypes.c_ulong),
                ("nFont", ctypes.c_ulong),
                ("dwFontSize", COORD),
                ("FontFamily", ctypes.c_uint),
                ("FontWeight", ctypes.c_uint),
                ("FaceName", ctypes.c_wchar * LF_FACESIZE)]


font_info = CONSOLE_FONT_INFOEX()
font_info.cbSize = ctypes.sizeof(CONSOLE_FONT_INFOEX)
font_info.nFont = 0
font_info.dwFontSize.X = 0
font_info.dwFontSize.Y = 22  # Increase font size
font_info.FontFamily = 54  # Set font family (use Monospace)

ctypes.windll.kernel32.SetCurrentConsoleFontEx(
    ctypes.windll.kernel32.GetStdHandle(STD_OUTPUT_HANDLE), False, ctypes.pointer(font_info)
)


np.set_printoptions(precision=15, suppress=True, linewidth=2000)

print("Enter the matrix you want to work with.")
main_matrix = request_matrix_interface("M")
print(f"The introduced matrix M is:\n{format_print_matrix(main_matrix)}")
print("-------------------------------------------------------------------------------")
e_1 = Matrix([[1, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]])

menu_options = {
    1: 'Matrix norm.',
    2: "Know if the matrix is a Mueller matrix.",
    3: 'Know if the matrix is K-irreducible.',
    4: 'Know if the matrix is K-primitive.',
    5: 'Approximation by an invertible matrix.',
    6: 'Approximation by a Mueller matrix.',
    7: 'Approximation by an invertible Mueller matrix.',
    8: 'Approximation by a K-primitive matrix.',
    9: 'Eigenvalue Calibration Method (ECM).',
    10: 'Enter a new matrix.',
    11: 'Exit.',
}


def print_menu():
    for key in menu_options.keys():
        print(key, '--', menu_options[key])


def option1():
    print("-----------------------------------")
    print(f'The norm of the introduced matrix M is: {matrix_norm(main_matrix)}')
    print("-----------------------------------")
    print('Press enter to return to the menu.')
    input()


def option2():
    plt.close("all")
    np.set_printoptions(precision=10, suppress=True, linewidth=2000)
    qm_0, mi_qm_0, minimum_qm_0, first_column_is_stokes, fig_0 = know_if_is_mueller(main_matrix)
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
    mat_1 = main_matrix * min_stock
    print(f"The result of the multiplication of M and v is:\n {np.array(mat_1).astype(np.float64)} ")

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
    print('Close the graph window and press enter to return to the menu.')
    input()


def option3():
    print("------------------------------------------------------------------------------------")
    k_irreducible("M", main_matrix)
    print("------------------------------------------------------------------------------------")
    print('Press enter to return to the menu.')
    input()


def option4():
    print("------------------------------------------------------------------------------------")
    k_primitive("M", main_matrix)
    print("------------------------------------------------------------------------------------")
    print('Press enter to return to the menu.')
    input()


def option5():
    print("------------------------------------------------------------------------------------")
    invertible_main_matrix, det_main_matrix = make_invertible("M", main_matrix)
    if det_main_matrix != 0:
        print(f"Then, M(inv)=M is: \n{np.array(invertible_main_matrix).astype(np.float64)}")
    else:
        print(f"The matrix M(inv) is: \n{np.array(invertible_main_matrix).astype(np.float64)}")
        print(f"The determinant of M(inv) is {invertible_main_matrix.det()}.")
    print("------------------------------------------------------------------------------------")
    print('Press enter to return to the menu.')
    input()


def option6():
    global main_matrix
    print("------------------------------------------------------------------------------------")
    print("------------------------------------------------------------------------------------")
    print(f"The input matrix M is:\n{np.array(main_matrix).astype(np.float64)}")
    print("------------------------------------------------------------------------------------")
    appr_mueller_matrix("M", main_matrix)
    print("------------------------------------------------------------------------------------")
    print('Press enter to return to the menu.')
    input()


def option7():
    print("------------------------------------------------------------------------------------")
    appr_invert_mueller_matrix("M", main_matrix)
    print("------------------------------------------------------------------------------------")
    print('Press enter to return to the menu.')
    input()


def option8():
    print("------------------------------------------------------------------------------------")
    appr_k_primitive("M", main_matrix)
    print("------------------------------------------------------------------------------------")
    print('Press enter to return to the menu.')
    input()


def option9():
    print("------------------------------------------------------------------------------------")
    try:
        m_appr = ecm(main_matrix)
        print(np.array(m_appr).astype(np.float64))
        print("This is the final matrix obtained by the ECM.")
    except Exception as e:
        print(f"An error occurred: {e}")
    print("------------------------------------------------------------------------------------")
    print('Press enter to return to the menu.')
    input()


def option10():
    global main_matrix
    main_matrix = request_matrix_interface("M")
    print(f"The input matrix M is:\n{format_print_matrix(main_matrix)}")


if __name__ == '__main__':
    while True:
        print_menu()
        option = ''
        try:
            option = int(input('Insert the number of the operation you want to perform in the matrix: '))
            print("------------------------------------------------------------------------------------")
        except:
            print('Invalid option. enter a number from 1 to 11.')
        # Check what choice was entered and act accordingly
        if option == 1:
            option1()
        elif option == 2:
            option2()
        elif option == 3:
            option3()
        elif option == 4:
            option4()
        elif option == 5:
            option5()
        elif option == 6:
            option6()
        elif option == 7:
            option7()
        elif option == 8:
            option8()
        elif option == 9:
            option9()
        elif option == 10:
            option10()
        elif option == 11:
            print("Thanks.")
            exit()
        else:
            print("-----------------------------------------------")
            print('Invalid option. enter a number from 1 to 11.')
            print("-----------------------------------------------")

