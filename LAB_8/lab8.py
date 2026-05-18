import math
import numpy as np

def f(x):
    return math.sin(x) - x / 2

def df(x):
    return math.cos(x) - 0.5

def d2f(x):
    return -math.sin(x)

def simple_iteration(x0, alpha, eps):
    x = x0
    iters = 0
    while True:
        x_new = x + alpha * f(x)
        iters += 1
        if abs(f(x_new)) < eps and abs(x_new - x) < eps:
            return x_new, iters
        x = x_new

def newton_method(x0, eps):
    x = x0
    iters = 0
    while True:
        x_new = x - f(x) / df(x)
        iters += 1
        if abs(f(x_new)) < eps and abs(x_new - x) < eps:
            return x_new, iters
        x = x_new

def chebyshev_method(x0, eps):
    x = x0
    iters = 0
    while True:
        fx = f(x)
        dfx = df(x)
        d2fx = d2f(x)
        x_new = x - fx / dfx - (d2fx / (2 * dfx)) * (fx / dfx)**2
        iters += 1
        if abs(f(x_new)) < eps and abs(x_new - x) < eps:
            return x_new, iters
        x = x_new

def secant_method(x0, x1, eps):
    iters = 0
    while True:
        fx1 = f(x1)
        fx0 = f(x0)
        x_new = x1 - fx1 * (x1 - x0) / (fx1 - fx0)
        iters += 1
        if abs(f(x_new)) < eps and abs(x_new - x1) < eps:
            return x_new, iters
        x0, x1 = x1, x_new

def parabola_method(x0, x1, x2, eps):
    iters = 0
    while True:
        f0, f1, f2 = f(x0), f(x1), f(x2)
        h1 = x1 - x0
        h2 = x2 - x1
        delta1 = (f1 - f0) / h1
        delta2 = (f2 - f1) / h2
        d = (delta2 - delta1) / (h2 + h1)
        b = delta2 + h2 * d
        D = math.sqrt(max(0, b**2 - 4 * f2 * d))
        
        if abs(b + D) > abs(b - D):
            E = b + D
        else:
            E = b - D
            
        h = -2 * f2 / E
        x_new = x2 + h
        iters += 1
        if abs(f(x_new)) < eps and abs(x_new - x2) < eps:
            return x_new, iters
        x0, x1, x2 = x1, x2, x_new

def inverse_interpolation(x0, x1, x2, eps):
    iters = 0
    while True:
        y0, y1, y2 = f(x0), f(x1), f(x2)
        x_new = (y1*y2)/((y0-y1)*(y0-y2))*x0 + \
                (y0*y2)/((y1-y0)*(y1-y2))*x1 + \
                (y0*y1)/((y2-y0)*(y2-y1))*x2
        iters += 1
        if abs(f(x_new)) < eps and abs(x_new - x2) < eps:
            return x_new, iters
        x0, x1, x2 = x1, x2, x_new

# Polynomial operations
def read_polynomial(filename):
    with open(filename, 'r') as file:
        coeffs = [float(x) for x in file.read().split()]
    return coeffs

def horner_eval(coeffs, x):
    val = coeffs[0]
    dval = coeffs[0]
    for a in coeffs[1:-1]:
        val = val * x + a
        dval = dval * x + val
    val = val * x + coeffs[-1]
    return val, dval

def newton_horner(coeffs, x0, eps):
    x = x0
    iters = 0
    while True:
        iters += 1
        val, dval = horner_eval(coeffs, x)
        if dval == 0:
            break
        x_new = x - val / dval
        if abs(val) < eps and abs(x_new - x) < eps:
            return x_new, iters
        x = x_new
        if iters > 1000:
            print("Newton Horner not converging")
            return x, iters
    return x, iters

def lin_method(coeffs, p0, q0, eps):
    n = len(coeffs) - 1
    p, q = p0, q0
    iters = 0
    while True:
        iters += 1
        c = [0] * (n + 1)
        c[0] = coeffs[0]
        c[1] = coeffs[1] - p * c[0]
        for k in range(2, n + 1):
            c[k] = coeffs[k] - p * c[k-1] - q * c[k-2]
        
        R = c[n-1]
        S = c[n]
        
        b = [0] * (n + 1)
        b[0] = c[0]
        b[1] = c[1] - p * b[0]
        for k in range(2, n):
            b[k] = c[k] - p * b[k-1] - q * b[k-2]
            
        det = b[n-2]**2 - b[n-3] * (b[n-1] - b[n-3]*q)
        if det == 0:
            print("Zero det in Bairstow")
            break
        dp = (c[n-1] * b[n-2] - c[n] * b[n-3]) / det
        dq = (c[n] * b[n-2] - c[n-1] * (b[n-1] - b[n-3]*q)) / det
        
        p += dp
        q += dq
        
        if abs(dp) < eps and abs(dq) < eps:
            break
        if iters > 1000:
            print("Bairstow not converging")
            break
            
    # Roots of x^2 + px + q = 0
    discriminant = p**2 - 4*q
    if discriminant < 0:
        real_part = -p / 2
        imag_part = math.sqrt(-discriminant) / 2
        return complex(real_part, imag_part), complex(real_part, -imag_part), iters
    else:
        r1 = (-p + math.sqrt(discriminant)) / 2
        r2 = (-p - math.sqrt(discriminant)) / 2
        return r1, r2, iters

def main():
    print("Лабораторна робота №8")
    a, b = -2, 2
    h = 0.1
    print("1. Табуляція функції f(x) = sin(x) - x/2:")
    
    # Point 1: Tabulate and find approximate roots
    nodes = []
    with open("tabulation.txt", "w") as f_out:
        for i in range(int((b - a) / h) + 1):
            x = a + i * h
            val = f(x)
            nodes.append((x, val))
            f_out.write(f"{x:.2f} {val:.6f}\n")
    print("Результати збережено в tabulation.txt\n")

    approx_roots = []
    for i in range(len(nodes) - 1):
        x_curr, y_curr = nodes[i]
        x_next, y_next = nodes[i+1]
        # Avoid double-counting if a node is exactly zero
        if y_curr * y_next < 0 or y_curr == 0:
            approx_roots.append((x_curr + x_next) / 2 if y_curr != 0 else x_curr)
            
    # Filter out duplicates that might appear if we land exactly on 0
    unique_roots = []
    for r in approx_roots:
        if not unique_roots or abs(r - unique_roots[-1]) > 0.1:
            unique_roots.append(r)
    approx_roots = unique_roots

    print(f"Наближені корені з табуляції: {approx_roots}")

    eps = 1e-5
    
    # First root (increasing, near 0)
    print("\n--- Корінь 1 (зростання, біля 0) ---")
    # approx_roots typically: ~ -1.89, 0.0, ~ 1.89
    # pick the one closest to 0
    x0_1 = min(approx_roots, key=lambda x: abs(x)) if approx_roots else 0.5
    print("Проста ітерація:", simple_iteration(x0_1, -1.0, eps))
    print("Ньютона:", newton_method(x0_1, eps))
    print("Чебишева:", chebyshev_method(x0_1, eps))
    print("Хорд:", secant_method(1.0, x0_1, eps))
    print("Парабол:", parabola_method(1.0, 0.8, x0_1, eps))
    print("Зворотня інтерполяція:", inverse_interpolation(1.0, 0.8, x0_1, eps))

    # Second root (decreasing, near 1.895)
    print("\n--- Корінь 2 (спадання, біля 1.895) ---")
    # pick the one closest to 1.895
    x0_2 = min(approx_roots, key=lambda x: abs(x - 1.895)) if approx_roots else 1.5
    print("Проста ітерація:", simple_iteration(x0_2, 1.0, eps))
    print("Ньютона:", newton_method(x0_2, eps))
    print("Чебишева:", chebyshev_method(x0_2, eps))
    print("Хорд:", secant_method(1.0, x0_2, eps))
    print("Парабол:", parabola_method(1.0, 1.2, x0_2, eps))
    print("Зворотня інтерполяція:", inverse_interpolation(1.0, 1.2, x0_2, eps))

    print("\n--- Алгебраїчне рівняння ---")
    coeffs = read_polynomial("coeffs.txt")
    print(f"Коефіцієнти: {coeffs}")
    
    # Point 5: Plotting the algebraic equation
    try:
        import matplotlib.pyplot as plt
        xs = np.linspace(-3, 3, 400)
        ys = [horner_eval(coeffs, x)[0] for x in xs]
        plt.figure()
        plt.plot(xs, ys, label="P(x)")
        plt.axhline(0, color='black', linewidth=1)
        plt.grid(True)
        plt.title("Графік алгебраїчного многочлена (Пункт 5)")
        plt.legend()
        plt.savefig("algebraic_plot.png")
        plt.close()
        print("Графік алгебраїчного многочлена збережено у 'algebraic_plot.png'")
    except ImportError:
        print("matplotlib не встановлено, пропускаємо побудову графіка.")

    root_real, iters_real = newton_horner(coeffs, 2.5, eps)
    print(f"Дійсний корінь (Ньютон+Горнер): {root_real:.6f}, ітерацій: {iters_real}")
    
    roots_complex = lin_method(coeffs, 0, 1, eps)
    print(f"Комплексні корені (метод Ліна/Баєрстоу): {roots_complex[0]}, {roots_complex[1]}, ітерацій: {roots_complex[2]}")

if __name__ == "__main__":
    main()
