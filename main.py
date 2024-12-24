from modules.method_z import *

if __name__ == "__main__":
    # Параметры задачи
    n, m = 3, 3  # размер сетки
    a, b, c, d = 0, 1, 0, 1  # границы области
    eps = 1e-6  # точность

    # Метод Зейделя
    max_iter = 10000
    v = method_z(max_iter, n, m, a, b, c, d, eps)

    # Вывод результата
    print("Решение:")
    print(v)
