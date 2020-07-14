#!/usr/bin/env python3
import numpy as np
import random

pi = 3.141592653


def rho_n(n):
    if n == 0:
        return n
    else:
        return n + omega/(pi * n)


def alpha_n(n):
    if n == 0:
        return pi/2.0
    else:
        return pi/2.0


def alpha_n0(n):
    if n == 0:
        return pi
    else:
        return pi/2.0


def set_bounds(size, step):
    array = []
    for i in range(size):
        array.append(i * step)
    return array


def f(x, t, eps):
    #  Подсчет суммы ряда при каждом фиксированном x и t, для заданного эпсилон
    n = 0
    gen = 1.0
    the_sum = 0.0
    while abs(gen) >= eps:
        gen = (((np.cos(rho_n(n) * x) * np.cos(rho_n(n) * t)) / alpha_n(n)) -
               ((np.cos(n * x) * np.cos(n * t)) / alpha_n0(n)))
        n += 1
        the_sum += gen
    return the_sum


def set_series_sum_matrix(bounds_array, eps):
    #  Заполнение матрицы ряда (значения ряда для всех х и t)
    matrix_series = []
    for i in range(len(bounds_array)):
        array = []
        for j in range(len(bounds_array)):
            x = bounds_array[i]
            t = bounds_array[j]
            array.append(f(x, t, eps))
        matrix_series.append(array)
    return matrix_series


def set_coefficients(bounds_array, series_array):
    #  задание коэффициентов СЛАУ, при каждом фиксированном x интегральное уравление дает некоторую СЛАУ.
    #  Проходимся по всем х и получаем N матриц и N столбцов свободных членов
    #  Функиця возвращает массив матриц и массив столбцов
    H = bounds_array[1] - bounds_array[0]
    matrix_array = []
    column_array = []
    #  заполнение вектора-столбца
    for i in range(len(bounds_array)):
        temp_col = []
        for j in range(i + 1):
            temp_col.append(-1.0 * series_array[i][j])
        column_array.append(temp_col)

    #  заполнение матрицы матриц
    for k in range(len(bounds_array)):
        temp_matrix = []
        temp_row = []
        if k == 0:
            temp_row.append(1.0)
            temp_matrix.append(temp_row)
        else:
            for i in range(k + 1):
                temp_row = []
                for j in range(k + 1):
                    if i == j:
                        if i == 0 or i == k:
                            temp_row.append(1.0 + series_array[j][i] * H / 2.0)
                        else:
                            temp_row.append(1.0 + series_array[j][i] * H)
                    else:
                        if j == 0 or j == k:
                            temp_row.append(series_array[j][i] * H / 2.0)
                        else:
                            temp_row.append(series_array[j][i] * H)
                temp_matrix.append(temp_row)
        matrix_array.append(temp_matrix)
    return matrix_array, column_array


def find_g(matrix_array, column_array):
    #  нахождение решения для каждой СЛАУ.
    i = 0
    solution_matrix = []
    for matrix in matrix_array:
        m = np.asarray(matrix)
        c = np.asarray(column_array[i])
        x = np.linalg.solve(m, c)
        i += 1
        x_list = x.tolist()
        solution_matrix.append(x_list)
    return solution_matrix


def find_q(g, step):
    array = []
    size = len(g[len(g) - 1])
    for i in range(size):
        if i == 0:
            temp = 2.0 * ((g[i + 1][i + 1] - g[i][i])/step)
        elif i == size - 1:
            temp = 2.0 * ((g[i][i] - g[i - 1][i - 1])/step)
        else:
            temp = 2.0 * ((g[i + 1][i + 1] - g[i - 1][i - 1])/(2.0 * step))
        array.append(temp)

    return array


def find_h(g):
    return g[0][0]


def find_H(h_low, q_array, step):
    integral = 0
    for i in range(len(q_array) - 1):
        temp = step * (q_array[i] + q_array[i + 1]) / 2.0
        integral += temp

    h_high = omega - h_low - 0.5 * integral

    return h_high


def exp_q(bounds_array):
    array = []
    for bound in bounds_array:
        array.append(2.0 * pow((2.0/pi - 1.0/pi), 2.0) / pow((1.0 + bound * (2.0/pi - 1.0/pi)), 2.0))

    return array


def exp_g(bounds_array):
    array = []
    for i in range(len(bounds_array)):
        temp = []
        for j in range(i + 1):
            temp.append(-(2.0/pi - 1.0/pi) / (1.0 + bounds_array[i] * (2.0/pi - 1.0/pi)))
        array.append(temp)

    return array


def exp_h():
    return -(2.0/pi - 1.0/pi)


def exp_H():
    return (2.0/pi - 1.0/pi) / (1.0 + (2.0/pi - 1.0/pi) * pi)


print("Enter number of bounds:")
N = int(input())

print("Enter epsilon:")
epsilon = float(input())
omega = 0.0
a = 0.0
b = pi
grid_step = (b - a) / N
N += 1
print("Grid step = ", grid_step)

bounds = set_bounds(N, grid_step)

file = open('output.txt', 'w')
file_q = open('q(x).txt', 'w')
for i in range(len(bounds)):
    file.write(str(i) + ")" + str(round(bounds[i], 5)) + '\n')

file.write('\n')

series_matrix = set_series_sum_matrix(bounds, epsilon)

file.write("F(x,t)\n")
for i in range(len(series_matrix)):
    for j in range(len(series_matrix[i])):
        file.write(str(round(series_matrix[i][j], 5)) + '\t')
    file.write('\n')

#  массивы матриц и столбцов свободных членов
matrix_vector, column_vector = set_coefficients(bounds, series_matrix)

file.write("\nМатрицы коэффициентов: \n")
counter = 0
for matrix_single in matrix_vector:
    sub_counter = 0
    file.write("\nMatrix[" + str(counter) + "]\n")
    for row in matrix_single:
        for i in range(len(row)):
            file.write(str(round(row[i], 5)) + " ")
        file.write(str(round(column_vector[counter][sub_counter], 5)) + '\n')
        sub_counter += 1
    counter += 1

g_solution = find_g(matrix_vector, column_vector)

file.write("\nФункция G(x,t): \n")
for row in g_solution:
    for i in range(len(row)):
        file.write(str(round(row[i], 6)) + " ")
    file.write('\n')

file.write("\nТочная функция G(x,t): \n")
for row in exp_g(bounds):
    for i in range(len(row)):
        file.write(str(round(row[i], 6)) + " ")
    file.write('\n')


q_solution = find_q(g_solution, grid_step)

file.write("\nФункция q(x): \n")
index = 0
for num in q_solution:
    file.write(str(round(num, 6)) + " ")
    file_q.write(str(round(bounds[index], 7)).replace('.', ',') + " " + str(round(num, 7)).replace('.', ',') + '\n')
    index += 1
file.write('\n')
file_q.close()

file.write("\nТочная функция q(x): \n")
for num in exp_q(bounds):
    file.write(str(round(num, 6)) + " ")
file.write('\n')

h = find_h(g_solution)
H = find_H(h, q_solution, grid_step)

file.write("\nПараметр h: \n" + str(round(h, 8)) + '\n')
file.write("Точный параметр h: \n" + str(round(exp_h(), 8)) + '\n')
file.write("\nПараметр H: \n" + str(round(H, 8)) + '\n')
file.write("Точный параметр H: \n" + str(round(exp_H(), 8)) + '\n')

file.write("Погрешности:\n")
file.write("Погрешность для h:" + str(round(abs(h - exp_h()), 4)) + '\n')
file.write("Погрешность для H:" + str(round(abs(H - exp_H()), 4)) + '\n')

delta = list()
for i in range(len(q_solution)):
    delta.append(abs(q_solution[i] - exp_q(bounds)[i]))

file.write("Наибольшая погрешность для функции q(x): " + str(round(max(delta), 10)) + '\n')

for i in range(len(g_solution)):
    for j in range(i):
        delta.append(abs(g_solution[i][j] - exp_g(bounds)[i][j]))

file.write("Наибольшая погрешность для функции G(x,t): " + str(round(max(delta), 4)) + '\n')

# ---------------ДОБАВЛЕНИЯ НА УСТОЙЧИВОСТЬ----------------------
for i in range(len(g_solution)):
    for j in range(len(g_solution[i])):
        z = (random.randint(0, 100) / 1000)
        g_solution[i][j] += z
        print(z)

file.write("\nФункция G(x,t): \n")
for row in g_solution:
    for i in range(len(row)):
        file.write(str(round(row[i], 6)) + " ")
    file.write('\n')

q_solution = find_q(g_solution, grid_step)

file.write("\nФункция q(x): \n")
for num in q_solution:
    file.write(str(round(num, 6)) + " ")
file.write('\n')

delta.clear()
for i in range(len(q_solution)):
    delta.append(abs(q_solution[i] - exp_q(bounds)[i]))

file.write("Наибольшая погрешность для функции q(x): " + str(round(max(delta), 10)) + '\n')
file.close()
