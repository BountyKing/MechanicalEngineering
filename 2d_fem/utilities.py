import numpy as np


def ksi(x, y, n1, n2, n3):
    x1, y1 = n1
    x2, y2 = n2
    x3, y3 = n3
    A = 0.5 * ((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1))
    a = (x3 * y1 - x1 * y3) / (2 * A)
    b = (y3 - y1) / (2 * A)
    c = (x1 - x3) / (2 * A)
    return a + b * x + c * y


def eta(x, y, n1, n2, n3):
    x1, y1 = n1
    x2, y2 = n2
    x3, y3 = n3
    A = 0.5 * ((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1))
    a = (x1 * y2 - x2 * y1) / (2 * A)
    b = (y1 - y2) / (2 * A)
    c = (x2 - x1) / (2 * A)
    return a + b * x + c * y


def N1(_ksi, _eta):
    return 1 - _ksi - _eta


def grad_N1():
    return np.array([-1, -1])


def N2(_ksi, _eta):
    return _ksi


def grad_N2():
    return np.array([1, 0])


def N3(_ksi, _eta):
    return _eta


def grad_N3():
    return np.array([0, 1])


def Jk(n1, n2, n3):
    x1, y1 = n1
    x2, y2 = n2
    x3, y3 = n3
    return np.array([[x2 - x1, x3 - x1],
                     [y2 - y1, y3 - y1]])


def Jk_1(n1, n2, n3):
    x1, y1 = n1
    x2, y2 = n2
    x3, y3 = n3
    A = 0.5 * ((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1))
    return np.array([[y3 - y1, x1 - x3],
                     [y1 - y2, x2 - x1]]) / (2 * A)


def get_Fe(nodes, func):
    f = np.empty(len(nodes))
    for i, n in enumerate(nodes):
        f[i] = func(n[0], n[1])
    return f


def B(n1, n2, n3):
    x1, y1 = n1
    x2, y2 = n2
    x3, y3 = n3
    A = 0.5 * ((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1))
    b1 = (y3 - y1) / (2 * A)
    c1 = (x3 - x2) / (2 * A)
    b2 = (y1 - y2) / (2 * A)
    c2 = (x1 - x3) / (2 * A)
    b3 = (y2 - y3) / (2 * A)
    c3 = (x2 - x1) / (2 * A)
    return (1 / (2 * A)) * np.array([[b1, b2, b3],
                                     [c1, c2, c3]])
