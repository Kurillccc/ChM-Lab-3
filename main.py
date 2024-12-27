import matplotlib.pyplot as plt
import tkinter as tk
import customtkinter as ctk

from tkinter import ttk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from tkinter import scrolledtext
from modules.method_z import *

def true_solution(x, y):
    return x**3 + y**2 + 3 # 1 - (x - 1)**2 - (y - 1/2)**2 # x**3 + y**2 + 3
# ------- Дельта от u -------
def f(xi, yj):
    return  - 6 * xi - 2 # -4

def mu1(yj, a):
    return true_solution(a, yj)

def mu2(yj, b):
    return true_solution(b, yj)

def mu3(xi, c):
    return true_solution(xi, c)

def mu4(xi, d):
    return true_solution(xi, d)

def buildMatrix(n, m, a, b, c, d):
    matrix = np.zeros(((n - 1) * (m - 1), (n - 1) * (m - 1)))
    vec = np.zeros((n - 1) * (m - 1))

    h = (b - a) / n
    k = (d - c) / m

    A = -2 * (1 / h ** 2 + 1 / k ** 2)
    for j in range(1, m):
        for i in range(1, n):
            numberEq = (j - 1) * (n - 1) + i - 1
            matrix[numberEq][numberEq] = A
            vec[numberEq] -= f(i * h, j * k)

            if i == 1:
                vec[numberEq] -= mu1(j * k, a) / h ** 2
            else:
                matrix[numberEq][numberEq - 1] = 1 / h ** 2

            if i == n - 1:
                vec[numberEq] -= mu2(j * k, b) / h ** 2
            else:
                matrix[numberEq][numberEq + 1] = 1 / h ** 2

            if j == 1:
                vec[numberEq] -= mu3(i * h, c) / k ** 2
            else:
                matrix[numberEq][numberEq - (n - 1)] = 1 / k ** 2

            if j == m - 1:
                vec[numberEq] -= mu4(i * h, d) / k ** 2
            else:
                matrix[numberEq][numberEq + (n - 1)] = 1 / k ** 2

    return matrix, vec

def show_matrix_vector(matrix, vec):
    matrix_vector_window = tk.Toplevel(root)
    matrix_vector_window.title("Матрица")

    matrix_frame = ctk.CTkFrame(matrix_vector_window)
    matrix_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
    vector_frame = ctk.CTkFrame(matrix_vector_window)
    vector_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

    matrix_label = ctk.CTkLabel(matrix_frame, text="Матрица")
    matrix_label.pack(pady=2)
    matrix_display = scrolledtext.ScrolledText(matrix_frame, width=70, height=10)
    matrix_display.pack(pady=5)
    matrix_text = '\n'.join(['\t'.join([f"{item:.2f}" for item in row]) for row in matrix])
    matrix_display.insert(tk.INSERT, matrix_text)

    vec_label = ctk.CTkLabel(vector_frame, text="Правая часть")
    vec_label.pack(pady=2)
    vec_display = scrolledtext.ScrolledText(vector_frame, width=10, height=10)
    vec_display.pack(pady=5)
    vec_text = '\n'.join([f"{item:.2f}" for item in vec])
    vec_display.insert(tk.INSERT, vec_text)

def show_table_window(solution, n, m, a, b, c, d):
    eps = float(eps_entry.get())
    Nmax = int(nmax_entry.get())

    table_window = tk.Toplevel(root)
    table_window.title("Таблица результатов")

    columns = ('i', 'xi', 'j', 'yj', 'Vij', 'Uij', 'diff')
    table = ttk.Treeview(table_window, columns=columns, show='headings')
    table.heading('i', text='i')
    table.heading('xi', text='xi')
    table.heading('j', text='j')
    table.heading('yj', text='yj')
    table.heading('Vij', text='Vij')
    table.heading('Uij', text='Uij')
    table.heading('diff', text='|Vij - Uij|')

    scrollbar = ttk.Scrollbar(table_window, orient=tk.VERTICAL, command=table.yview)
    table.configure(yscrollcommand=scrollbar.set)
    scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

    h = (b - a) / n
    k = (d - c) / m
    max_diff = 0

    difff = []

    for j in range(m + 1):
        yj = c + j * k
        for i in range(n + 1):
            xi = a + i * h
            Vij = solution[i, j]
            Uij = true_solution(xi, yj)
            diff = abs(Vij - Uij)
            difff.append(diff)
            max_diff = max(max_diff, diff)
            table.insert('', tk.END, values=(i, "{:.5f}".format(xi), j, "{:.5f}".format(yj), "{:.20f}".format(Vij), "{:.20f}".format(Uij), "{:.20f}".format(diff)))

    matrix, vec = buildMatrix(n, m, a, b, c, d)
    result, x_for_pogr = seidel_method(matrix, vec, eps, Nmax)
    number = len(x_for_pogr)

    table.pack(expand=True, fill='both')

    max_diff_label = ctk.CTkLabel(table_window,
                                  text=f"Норма погрешность: {max_diff}\
                                  \n ||Z^(N)||2 = {np.sqrt(np.dot(difff,difff))}\
                                  \nТочность на выходе (EpsN) = {max(x_for_pogr[number - 1] - x_for_pogr[number - 2])}\
                                  \nНорма невязка на выходе ||R^(N)||2 = {np.sqrt(np.dot(np.dot(matrix, result) - vec, np.dot(matrix, result) - vec))}",
                                  text_color='black')
    max_diff_label.pack()

# ---------------------- Функции для кнопок ----------------------
def update():
    a = int(a_entry.get())
    b = int(b_entry.get())
    c = int(c_entry.get())
    d = int(d_entry.get())
    n = int(n_entry.get())
    m = int(m_entry.get())
    eps = float(eps_entry.get())
    Nmax = int(nmax_entry.get())

    matrix, vec = buildMatrix(n, m, a, b, c, d)
    result, x_for_pogr = seidel_method(matrix, vec, eps, Nmax)


    solution = np.zeros((n + 1, m + 1))
    solution[:, 0] = [mu3(xi, c) for xi in np.linspace(a, b, n + 1)]
    solution[:, -1] = [mu4(xi, d) for xi in np.linspace(a, b, n + 1)]
    solution[0, :] = [mu1(yj, a) for yj in np.linspace(c, d, m + 1)]
    solution[-1, :] = [mu2(yj, b) for yj in np.linspace(c, d, m + 1)]

    for j in range(1, m):
        for i in range(1, n):
            solution[i, j] = result[(j - 1) * (n - 1) + (i - 1)]

    x = np.linspace(a, b, n + 1)
    y = np.linspace(c, d, m + 1)
    X, Y = np.meshgrid(x, y)

    ax.clear()
    ax.plot_surface(X, Y, solution.T, cmap='viridis')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('U(x,y)')
    plt.title('Численное решение')

    canvas.draw()

def results_table_show():
    a = int(a_entry.get())
    b = int(b_entry.get())
    c = int(c_entry.get())
    d = int(d_entry.get())
    n = int(n_entry.get())
    m = int(m_entry.get())
    eps = float(eps_entry.get())
    Nmax = int(nmax_entry.get())

    matrix, vec = buildMatrix(n, m, a, b, c, d)
    result, x_for_pogr = seidel_method(matrix, vec, eps, Nmax)

    solution = np.zeros((n + 1, m + 1))
    solution[:, 0] = [mu3(xi, c) for xi in np.linspace(a, b, n + 1)]
    solution[:, -1] = [mu4(xi, d) for xi in np.linspace(a, b, n + 1)]
    solution[0, :] = [mu1(yj, a) for yj in np.linspace(c, d, m + 1)]
    solution[-1, :] = [mu2(yj, b) for yj in np.linspace(c, d, m + 1)]

    for j in range(1, m):
        for i in range(1, n):
            solution[i, j] = result[(j - 1) * (n - 1) + (i - 1)]
    show_table_window(solution, n, m, a, b, c, d)


def results_table_show_modern(index):
    a = int(a_entry.get())
    b = int(b_entry.get())
    c = int(c_entry.get())
    d = int(d_entry.get())
    n = int(n_entry.get())
    m = int(m_entry.get())
    eps = float(eps_entry.get())
    Nmax = int(nmax_entry.get())

    # Индекс уменьшаем на 1, так как индексы могут начинаться с 0
    index = index - 1

    matrix, vec = buildMatrix(n, m, a, b, c, d)
    x_new, x_layer, s = seidel_method_modern(matrix, vec, eps, Nmax, index)

    plt.figure(figsize=(8, 5))

    # Если индекс равен -1, рисуем все слои
    if index == -1:
        for i, layer in enumerate(x_layer):
            plt.plot(range(1, len(layer) + 1), layer, marker='o', label=f"Слой {i + 1}")
    else:
        # Рисуем только один слой
        plt.plot(range(1, len(x_new) + 1), x_new, marker='o', label=f"{index} слой")

    plt.title(f"Сходимость достигнута за {s} итераций.")
    plt.xlabel("Индекс переменной")
    plt.ylabel("X")
    plt.grid(True)
    plt.legend()

    plt.show()

def matrix_show():
    a = int(a_entry.get())
    b = int(b_entry.get())
    c = int(c_entry.get())
    d = int(d_entry.get())
    n = int(n_entry.get())
    m = int(m_entry.get())

    matrix, vec = buildMatrix(n, m, a, b, c, d)
    show_matrix_vector(matrix, vec)

def Layer():
    plt.ion()
    index = int(index_entry.get())
    results_table_show_modern(index + 1)


if __name__ == "__main__":
    root = tk.Tk()
    root.title("ЧМ лабораторная работа 3")

    frame = ttk.Frame(root)
    frame.pack(padx=5, pady=5)

    root.geometry("1100x700")
    root.minsize(1100, 700)

    window_width = 1100
    window_height = 900

    screen_width = root.winfo_screenwidth()
    screen_height = root.winfo_screenheight()

    position_top = int(screen_height / 2 - window_height / 2)
    position_right = int(screen_width / 2 - window_width / 2)

    root.geometry(f'{window_width}x{window_height}+{position_right}+{position_top}')

    # ---------------------- Входные параметры ----------------------
    a_label = ttk.Label(frame, text="a")
    a_label.grid(row=0, column=0)
    a_entry = ttk.Entry(frame)
    a_entry.grid(row=0, column=1)
    a_entry.insert(0, "0")

    b_label = ttk.Label(frame, text="b")
    b_label.grid(row=1, column=0)
    b_entry = ttk.Entry(frame)
    b_entry.grid(row=1, column=1)
    b_entry.insert(0, "1")

    c_label = ttk.Label(frame, text="c")
    c_label.grid(row=2, column=0)
    c_entry = ttk.Entry(frame)
    c_entry.grid(row=2, column=1)
    c_entry.insert(0, "0")

    d_label = ttk.Label(frame, text="d")
    d_label.grid(row=3, column=0)
    d_entry = ttk.Entry(frame)
    d_entry.grid(row=3, column=1)
    d_entry.insert(0, "1")

    n_label = ttk.Label(frame, text="n")
    n_label.grid(row=0, column=2)
    n_entry = ttk.Entry(frame)
    n_entry.grid(row=0, column=3)
    n_entry.insert(0, "3")

    m_label = ttk.Label(frame, text="m")
    m_label.grid(row=1, column=2)
    m_entry = ttk.Entry(frame)
    m_entry.grid(row=1, column=3)
    m_entry.insert(0, "3")

    eps_label = ttk.Label(frame, text="Eps.")
    eps_label.grid(row=0, column=4)
    eps_entry = ttk.Entry(frame)
    eps_entry.grid(row=0, column=5)
    eps_entry.insert(0, "1e-8")

    nmax_label = ttk.Label(frame, text="Nmax")
    nmax_label.grid(row=1, column=4)
    nmax_entry = ttk.Entry(frame)
    nmax_entry.grid(row=1, column=5)
    nmax_entry.insert(0, "10000")

    index_label = ttk.Label(frame, text="     index \n('-1' for all)")
    index_label.grid(row=0, column=6)
    index_entry = ttk.Entry(frame)
    index_entry.grid(row=0, column=7)
    index_entry.insert(0, "1")

    # ---------------------- Кнопки ----------------------
    update_button = ttk.Button(frame, text="Обновить", command=update)
    update_button.grid(row=11, column=1)

    update_button = ttk.Button(frame, text="Таблица результатов", command=results_table_show)
    update_button.grid(row=11, column=3)

    update_button = ttk.Button(frame, text="Матрица", command=matrix_show)
    update_button.grid(row=11, column=5)

    update_button = ttk.Button(frame, text="Слой", command=Layer)
    update_button.grid(row=1, column=7)

    entries = {}

    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    canvas = FigureCanvasTkAgg(fig, master=root)
    canvas_widget = canvas.get_tk_widget()
    canvas_widget.pack(fill=tk.BOTH, expand=True)

    update()

    root.mainloop()