import numpy as np


# 1d Test functions
def func_1d_1(x):
    return x ** 2


def func_1d_2(x):
    return x ** 3 / 10**3


def func_1d_3(x):
    return x ** 4 / 10**4


def func_1d_4(x):
    return x ** 5 / 10**5


def func_1d_5(x):
    return np.sin(x)


def func_1d_6(x):
    return np.exp(x)


func_list_1d = [func_1d_1, func_1d_2, func_1d_3, func_1d_4, func_1d_5, func_1d_6]
func_str_list_1d = ["f(x) = x^2", "f(x) = x^3", "f(x) = x^4", "f(x) = x^5", "f(x) = sin(x)", "f(x) = exp(x)"]
