from scipy import integrate
import numpy
from plot import draw_plot

rho = 1
a = 0
b = 3
length = b - a


def eps(x):
    if 0 <= x <= 1:
        return 10
    elif 1 < x <= 2:
        return 5
    else:
        return 1


def e(x, i, h):
    x0 = (i-1)*h
    x1 = (i+1)*h
    xi = i*h

    if x0 < x < xi:
        return (x-x0)/h
    elif xi <= x < x1:
        return (x1-x)/h
    else:
        return 0


def e_prim(x, i, h):
    x0 = (i - 1) * h
    x1 = (i + 1) * h
    xi = i * h

    if x0 < x < xi:
        return 1/h
    elif xi <= x < x1:
        return (-1)/h
    else:
        return 0


def B(i, j, h):
    part = e(i, 0, h) * e(j, 0, h)

    integral_upper_bound = min((i+1)*h, (j+1)*h, b)
    integral_lower_bound = max((i-1)*h, (j-1)*h, a)

    integral = integrate.quad(lambda x: e_prim(x, i, h) * e_prim(x, j, h), integral_lower_bound, integral_upper_bound)[0]
    return part - integral


def L(i, h):
    part = 5 * e(i, 0, h)

    integral_upper_bound = min((i+1)*h, b)
    integral_lower_bound = max((i-1)*h, a)

    integral = integrate.quad(lambda x: rho/eps(x)*e(x, i, h), integral_lower_bound, integral_upper_bound)[0]
    return part - integral


def l2(i, n, h):
    return L(i, h) - 2*B(n, i, h)


def create_matrix_b(n, h):
    matrix = []
    for i in range(n):
        matrix_row = []
        for j in range(n):
            if abs(i-j) > 1:
                matrix_row.append(0)
            else:
                matrix_row.append(B(i, j, h))
        matrix.append(matrix_row)
    return matrix


def create_matrix_l(n, h):
    matrix = []
    for i in range(n):
        matrix.append(l2(i, n, h))
    return matrix


def fem(n):
    numpy.set_printoptions(precision=4, floatmode='fixed')
    h = length / n

    matrix_b = numpy.array(create_matrix_b(n, h))
    matrix_l = numpy.array(create_matrix_l(n, h))

    result_x = [i*h for i in range(n+1)]
    result_matrix_y = numpy.linalg.solve(matrix_b, matrix_l)
    result = numpy.append(result_matrix_y, 2)

    # print(result_x)
    # print(result)

    draw_plot(result_x, result)
