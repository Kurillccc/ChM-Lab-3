import numpy as np

def seidel_method(A, b, eps, Nmax):
    count = len(b)
    x = np.zeros(count)
    for S in range(Nmax):
        x_new = np.copy(x)
        for i in range(count):
            s1 = sum(A[i][j] * x_new[j] for j in range(i))
            s2 = sum(A[i][j] * x[j] for j in range(i + 1, count))
            x_new[i] = (b[i] - s1 - s2) / A[i][i]

        if np.linalg.norm(x_new - x, ord=np.inf) < eps:
            return x_new

        x = x_new

    return x
# ------------- Для вывода слоя -------------
def seidel_method_modern(A, b, eps, Nmax, layer_index=None):
    count = len(b)
    x = np.zeros(count)  # Начальное приближение
    x_all_layers = []  # Список для всех слоев, если index == -1

    for S in range(Nmax):
        x_new = np.copy(x)
        for i in range(count):
            # Разбиваем сумму на две части: до текущего элемента и после
            s1 = sum(A[i][j] * x_new[j] for j in range(i))
            s2 = sum(A[i][j] * x[j] for j in range(i + 1, count))
            x_new[i] = (b[i] - s1 - s2) / A[i][i]

        # Если индекс задан, сохраняем слой
        if layer_index is not None and S == layer_index:
            x_layer = x_new
            print(f"Состояние на итерации {S}: {x_new}")
            return x_new, x_layer, S + 1

        # Сохраняем все значения x на каждой итерации, если layer_index == -1
        if layer_index == -1:
            x_all_layers.append(x_new.copy())
            print(f"Состояние на итерации {S}: {x_new}")

        # Проверка сходимости
        if np.linalg.norm(x_new - x, ord=np.inf) < eps:
            print(f"Сходимость достигнута за {S + 1} итераций.")
            return x_new, x_all_layers, S + 1

        x = x_new

    print("Сходимость не достигнута за максимальное количество итераций.")
    return x, x_all_layers