import numpy as np

# Функции
def f(x, y):
    return x**3 + y**2 + 2

# Граничные условия
def boundary_condition(x, y):
    return 0  # можно заменить на конкретные условия

def method_z(max_iter, n, m, a, b, c, d, eps):
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

    for iteration in range(max_iter):
        max_error = 0
        for i in range(1, n - 1):
            for j in range(1, m - 1):
                v_new = (v[i - 1, j] + v[i + 1, j]) / hx**2 + (v[i, j - 1] + v[i, j + 1]) / hy**2 - f(x[i], y[j])
                v_new /= 2 * (1 / hx**2 + 1 / hy**2)
                max_error = max(max_error, abs(v[i, j] - v_new))
                v[i, j] = v_new
        if max_error < eps:
            print(f"Сходимость достигнута за {iteration} итераций")
            return v
            break
    else:
        return v
        print("Сходимость не достигнута")