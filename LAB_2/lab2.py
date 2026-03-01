import csv
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

def read_data(filename):
    """Зчитує дані (X та Y) з CSV файлу."""
    x = []
    y = []
    with open(filename, 'r', newline='') as file:
        reader = csv.DictReader(file)
        fieldnames = reader.fieldnames
        for row in reader:
            x.append(float(row[fieldnames[0]]))
            y.append(float(row[fieldnames[1]]))
    return np.array(x), np.array(y)
# -------------------------------------------------------------------
def divided_diff_table(x, y):
    """"""
    n = len(y)
    coef = np.zeros([n, n])
    coef[:,0] = y
    for j in range(1, n):
        for i in range(n - j):
            coef[i][j] = (coef[i+1][j-1] - coef[i][j-1]) / (x[i+j] - x[i])
    return coef

def newton_poly_div(coef, x_data, x):
    """"""
    n = len(x_data) - 1
    p = coef[0][n]
    for k in range(1, n + 1):
        p = coef[0][n-k] + (x - x_data[n-k])*p
    return p

# -------------------------------------------------------------------
def finite_diff_table(y):
    """"""
    n = len(y)
    diffs = np.zeros([n, n])
    diffs[:,0] = y
    for j in range(1, n):
        for i in range(n - j):
            diffs[i][j] = diffs[i+1][j-1] - diffs[i][j-1]
    return diffs

def factorial_poly(diffs, y_data, q):
    """"""
    n = len(y_data)
    p = diffs[0][0]
    q_term = 1.0
    for k in range(1, n):
        q_term *= (q - k + 1) / k
        p += diffs[0][k] * q_term
    return p

# -------------------------------------------------------------------
def lagrange_poly(x_data, y_data, x):
    """"""
    total = 0
    n = len(x_data)
    for i in range(n):
        term = y_data[i]
        for j in range(n):
            if i != j:
                term = term * (x - x_data[j]) / (x_data[i] - x_data[j])
        total += term
    return total

# -------------------------------------------------------------------
def test_function_1(x):
    return np.sin(x)

def runge_function(x):
    return 1.0 / (1.0 + 25.0 * x**2)

def research_part():
    print("\n--- Дослідницька частина ---")
    
    print("1. Дослідження впливу кроку (інтервал [0, 6], f(x) = sin(x))")
    interval_1 = [0, 6]
    nodes_counts_1 = [5, 10, 20]
    linestyles = ['-', '--', ':']
    
    plt.figure(figsize=(10, 6))
    x_exact = np.linspace(interval_1[0], interval_1[1], 500)
    plt.plot(x_exact, test_function_1(x_exact), 'k-', label="Exact sin(x)", linewidth=2, alpha=0.5)
    
    for idx, n in enumerate(nodes_counts_1):
        x_nodes = np.linspace(interval_1[0], interval_1[1], n)
        y_nodes = test_function_1(x_nodes)
        coef = divided_diff_table(x_nodes, y_nodes)
        y_interp = [newton_poly_div(coef, x_nodes, xi) for xi in x_exact]
        error = np.max(np.abs(y_interp - test_function_1(x_exact)))
        print(f"   Вузлів: {n}, Макс. похибка: {error:.2e}")
        # Make the lines thicker for smaller n so the exact line is still visible under it, or vice versa
        plt.plot(x_exact, y_interp, label=f"Newton n={n}", linewidth=max(2, 5 - idx), linestyle=linestyles[idx], alpha=0.8)
    
    plt.title("Вплив кроку (фіксований інтервал [0, 6])")
    plt.legend()
    plt.grid()
    plt.savefig("step_inf.png")
    plt.close()
    
    print("2. Вплив кількості вузлів (крок h=0.5, f(x) = sin(x))")
    step = 0.5
    intervals_2 = [[-1, 1], [-2, 2], [-4, 4]]
    
    for interval in intervals_2:
        n = int((interval[1] - interval[0]) / step) + 1
        x_nodes = np.linspace(interval[0], interval[1], n)
        y_nodes = test_function_1(x_nodes)
        coef = divided_diff_table(x_nodes, y_nodes)
        
        x_test = np.linspace(interval[0], interval[1], 300)
        y_interp = [newton_poly_div(coef, x_nodes, xi) for xi in x_test]
        error = np.max(np.abs(y_interp - test_function_1(x_test)))
        print(f"   Інтервал: {interval}, Вузлів: {n}, Макс. похибка: {error:.2e}")

    print("3. Аналіз ефекту Рунге (f(x) = 1/(1+25x^2) на [-1, 1])")
    plt.figure(figsize=(10, 6))
    x_exact_runge = np.linspace(-1, 1, 500)
    plt.plot(x_exact_runge, runge_function(x_exact_runge), 'k--', label="1/(1+25x^2)", linewidth=5)
    
    for idx, n in enumerate([5, 10, 20]):
        x_nodes = np.linspace(-1, 1, n)
        y_nodes = runge_function(x_nodes)
        coef = divided_diff_table(x_nodes, y_nodes)
        y_interp = [newton_poly_div(coef, x_nodes, xi) for xi in x_exact_runge]
        plt.plot(x_exact_runge, y_interp, label=f"n={n}", linewidth=max(1, 3 - idx))
        
    plt.title("Аналіз ефекту Рунге")
    plt.ylim(-1, 2)
    plt.legend()
    plt.grid()
    plt.savefig("runge.png")
    plt.close()
    
    print("4. Порівняння Ньютона та Лагранжа на f(x)=sin(x)")
    x_nodes = np.linspace(0, 3, 5)
    y_nodes = test_function_1(x_nodes)
    coef = divided_diff_table(x_nodes, y_nodes)
    
    test_pt = 1.5
    val_newton = newton_poly_div(coef, x_nodes, test_pt)
    val_lagrange = lagrange_poly(x_nodes, y_nodes, test_pt)
    print(f"   Т. {test_pt}: Ньютон={val_newton:.6f}, Лагранж={val_lagrange:.6f}, Точне={np.sin(test_pt):.6f}")

# -------------------------------------------------------------------
def run_variant_5():
    csv_file = "data_var5.csv"
    predict_x = 1000
    
    with open(csv_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Objects', 'FPS'])
        writer.writerows([[100, 120], [200, 110], [400, 90], [800, 65], [1600, 40]])
    print(f"[+] Сгенеровано файл {csv_file}")

    x_data, y_data = read_data(csv_file)
    print(f"\n--- ВАРІАНТ 5. Аналіз FPS від Objects ---")
    print(f"X (Об'єкти): {x_data}")
    print(f"Y (FPS): {y_data}")
    
    coef_div = divided_diff_table(x_data, y_data)
    y_pred_newton = newton_poly_div(coef_div, x_data, predict_x)
    print(f"\n[Многочлен Ньютона] Прогноз для x = {predict_x}: {y_pred_newton:.4f} FPS")
    
    print("\nТаблиця розділених різниць (перші порядки):")
    for i in range(len(y_data)):
        row_str = f"x={x_data[i]:.0f}, y={y_data[i]:.4f} | "
        for j in range(1, len(y_data) - i):
            row_str += f"d^{j}y = {coef_div[i][j]:.8e} | "
        print(row_str)
        
    i_data = np.array([math.log2(x / 100) for x in x_data])
    diffs = finite_diff_table(y_data)
    
    q_pred = math.log2(predict_x / 100)
    y_pred_fact = factorial_poly(diffs, y_data, q=q_pred)
    print(f"\n[Факторіальні многочлени] Прогноз для x = {predict_x}: {y_pred_fact:.4f} FPS")
    
    max_objects = 100
    for test_obj in range(100, 2000):
        fps_pred = newton_poly_div(coef_div, x_data, test_obj)
        if fps_pred >= 60:
            max_objects = test_obj
        else:
            break
            
    print(f"\n[Аналіз] Кількість об'єктів для підтримки стабільних 60+ FPS повинна бути <= {max_objects}")

    plt.figure(figsize=(10, 6))
    plt.scatter(x_data, y_data, color='red', label='Вузли', zorder=5, s=60)
    
    x_dense = np.linspace(100, 1800, 200)
    y_interp = [newton_poly_div(coef_div, x_data, xi) for xi in x_dense]
    
    plt.plot(x_dense, y_interp, 'b-', label='Інтерполяція', linewidth=3)
    plt.scatter([predict_x], [y_pred_newton], color='green', s=120, label=f'Прогноз у т. {predict_x} ({y_pred_newton:.1f} FPS)', zorder=6)
    
    plt.axhline(60, color='gray', linestyle='--', label='60 FPS', linewidth=2)
    plt.axvline(max_objects, color='orange', linestyle='--', label=f'Ліміт: {max_objects} об\'єктів', linewidth=2)
    
    plt.title('Варіант 5: Прогноз FPS в процесі оптимізації рушія')
    plt.xlabel('Objects (Кількість об\'єктів)')
    plt.ylabel('FPS')
    plt.legend()
    plt.grid(True)
    plt.savefig("var5.png")
    plt.close()
    print("[+] Графік інтерполяції збережено як 'var5.png'")

if __name__ == "__main__":
    run_variant_5()
    research_part()
    print("\nВиконання лабораторної роботи успішно завершено.")
