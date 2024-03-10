import numpy as np


# 1d Test functions
def func_2d_1(x, y):
    return (x + y) ** 2


def func_3d_1(x, y, z):
    return (x + y + z) ** 2


def func_4d_1(x, y, z, q):
    return (x + y + z + q) ** 2


def func_5d_1(x, y, z, q, r):
    return (x + y + z + q + r) ** 2


def func_2d_2(x, y):
    return np.sin(x + y)


def func_2d_3(x, y):
    return np.exp(x + y)


func_list_nd = [(func_2d_1, 2), (func_2d_2, 2), (func_2d_3, 2), (func_3d_1, 3), (func_4d_1, 4), (func_5d_1, 5)]
func_desc_nd = [
    'f(x, y) = (x + y)^2',
    'f(x, y) = sin(x + y)',
    'f(x, y) = exp(x + y)',
    'f(x, y, z) = (x + y + z)^2',
    'f(x, y, z, q) = (x + y + z + q)^2',
    'f(x, y, z, q, r) = (x + y + z + q + r)^2'
]

