from sympy import *
from sympy import Rational
import numpy as np
import PySimpleGUI as psg


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


def request_matrix_interface(matrix_name):
    final_matrix = []
    l1 = psg.Text(f"Enter the values of the matrix {matrix_name}: ", key='-OUT-', font=('Arial Bold',
                                                                                        12), size=(80, 1),
                  justification='center')
    l11 = psg.Text('[1,1]')
    c11 = psg.Input('', enable_events=True, key='-INPUT-11',
                    font=('Arial Bold', 12), size=(12, 1), justification='left')
    l12 = psg.Text('[1,2]')
    c12 = psg.Input('', enable_events=True, key='-INPUT-12',
                    font=('Arial Bold', 12), size=(12, 1), justification='left')
    l13 = psg.Text('[1,3]')
    c13 = psg.Input('', enable_events=True, key='-INPUT-13',
                    font=('Arial Bold', 12), size=(12, 1), justification='left')
    l14 = psg.Text('[1,4]')
    c14 = psg.Input('', enable_events=True, key='-INPUT-14',
                    font=('Arial Bold', 12), size=(12, 1), justification='left')
    ###############################
    l21 = psg.Text('[2,1]')
    c21 = psg.Input('', enable_events=True, key='-INPUT-21',
                    font=('Arial Bold', 12), size=(12, 1), justification='left')
    l22 = psg.Text('[2,2]')
    c22 = psg.Input('', enable_events=True, key='-INPUT-22',
                    font=('Arial Bold', 12), size=(12, 1), justification='left')
    l23 = psg.Text('[2,3]')
    c23 = psg.Input('', enable_events=True, key='-INPUT-23',
                    font=('Arial Bold', 12), size=(12, 1), justification='left')
    l24 = psg.Text('[2,4]')
    c24 = psg.Input('', enable_events=True, key='-INPUT-24',
                    font=('Arial Bold', 12), size=(12, 1), justification='left')
    ###############################
    l31 = psg.Text('[3,1]')
    c31 = psg.Input('', enable_events=True, key='-INPUT-31',
                    font=('Arial Bold', 12), size=(12, 1), justification='left')
    l32 = psg.Text('[3,2]')
    c32 = psg.Input('', enable_events=True, key='-INPUT-32',
                    font=('Arial Bold', 12), size=(12, 1), justification='left')
    l33 = psg.Text('[3,3]')
    c33 = psg.Input('', enable_events=True, key='-INPUT-33',
                    font=('Arial Bold', 12), size=(12, 1), justification='left')
    l34 = psg.Text('[3,4]')
    c34 = psg.Input('', enable_events=True, key='-INPUT-34',
                    font=('Arial Bold', 12), size=(12, 1), justification='left')
    ###############################
    l41 = psg.Text('[4,1]')
    c41 = psg.Input('', enable_events=True, key='-INPUT-41',
                    font=('Arial Bold', 12), size=(12, 1), justification='left')
    l42 = psg.Text('[4,2]')
    c42 = psg.Input('', enable_events=True, key='-INPUT-42',
                    font=('Arial Bold', 12), size=(12, 1), justification='left')
    l43 = psg.Text('[4,3]')
    c43 = psg.Input('', enable_events=True, key='-INPUT-43',
                    font=('Arial Bold', 12), size=(12, 1), justification='left')
    l44 = psg.Text('[4,4]')
    c44 = psg.Input('', enable_events=True, key='-INPUT-44',
                    font=('Arial Bold', 12), size=(12, 1), justification='left')
    b1 = psg.Button('Accept', key='-OK-', font=('Arial Bold', 12))
    b2 = psg.Button('Exit', font=('Arial Bold', 12))
    layout = [[l1], [l11, c11, l12, c12, l13, c13, l14, c14],
              [l21, c21, l22, c22, l23, c23, l24, c24],
              [l31, c31, l32, c32, l33, c33, l34, c34],
              [l41, c41, l42, c42, l43, c43, l44, c44], [b1, b2]]
    window = psg.Window('Enter matrix', layout, size=(690, 200))
    while True:
        event, values = window.read()
        if event.startswith('-INPUT-'):
            if len(values[event]) > 0 and values[event][-1] not in '0123456789.-':
                psg.popup("Invalid input. Please enter only decimal numbers.")
                window[event].update(values[event][:-1])
        if event == '-OK-':
            final_matrix.clear()
            for item in values:
                if values[item] == '':
                    psg.popup("There is an empty entrance")
                    print(item + " empty")
                    break
                else:
                    try:
                        final_matrix.append(Rational(str(float(values[item]))))
                        if values['-INPUT-11'] == '0':
                            psg.popup("The value at coordinate (1, 1) must be non-zero. Please try again.")
                            break
                    except:
                        psg.popup("Invalid input. Please enter only decimal numbers")

            if len(final_matrix) == 16:
                final_matrix_np = np.array(final_matrix).reshape((4, 4))
                window.close()
                return Matrix(final_matrix_np)

        if event == psg.WIN_CLOSED or event == 'Exit':
            break
    exit()

