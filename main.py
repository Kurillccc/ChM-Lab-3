import numpy as np
from modules.method_z import *

# Функции
def f(x, y):
    return x**3 + y**2 + 3

# Граничные условия
def boundary_condition(x, y):
    return 0  # можно заменить на конкретные условия

if __name__ == "__main__":
    # Параметры задачи
    n, m = 3, 3  # размер сетки
    a, b, c, d = 0, 1, 0, 1  # границы области
    eps = 1e-6  # точность

    # Шаг сетки
    hx = (b - a) / (n - 1)
    hy = (d - c) / (m - 1)

    # Инициализация сетки
    x = np.linspace(a, b, n)
    y = np.linspace(c, d, m)
    v = np.zeros((n, m))  # начальное приближение

    # Применяем граничные условия
    for i in range(n):
        v[i, 0] = boundary_condition(x[i], c)
        v[i, -1] = boundary_condition(x[i], d)
    for j in range(m):
        v[0, j] = boundary_condition(a, y[j])
        v[-1, j] = boundary_condition(b, y[j])

    # Метод Зейделя
    max_iter = 10000
    method_z(max_iter)

    # Вывод результата
    print("Решение:")
    print(v)
