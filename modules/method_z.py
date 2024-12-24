def method_z(max_iter):
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
            break
    else:
        print("Сходимость не достигнута")