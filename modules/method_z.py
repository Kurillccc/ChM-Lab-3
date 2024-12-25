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