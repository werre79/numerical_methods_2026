import csv
import math
import matplotlib.pyplot as plt

def read_data(filename):
    x = []
    y = []
    with open(filename, 'r') as f:
        reader = csv.reader(f)
        next(reader) # skip header
        for row in reader:
            x.append(float(row[0]))
            y.append(float(row[1]))
    return x, y

def form_matrix(x, m):
    n = len(x)
    A = [[0.0] * (m + 1) for _ in range(m + 1)]
    for i in range(m + 1):
        for j in range(m + 1):
            s = 0.0
            for k in range(n):
                s += x[k] ** (i + j)
            A[i][j] = s
    return A

def form_vector(x, y, m):
    n = len(x)
    b = [0.0] * (m + 1)
    for i in range(m + 1):
        s = 0.0
        for k in range(n):
            s += y[k] * (x[k] ** i)
        b[i] = s
    return b

def gauss_solve(A, b):
    # n is actually m+1 here
    n = len(b)
    
    # Create copies to not mutate original
    A = [row[:] for row in A]
    b = b[:]
    
    # Forward pass with partial pivoting
    for k in range(n - 1):
        # find max element in column k
        max_val = abs(A[k][k])
        max_row = k
        for i in range(k + 1, n):
            if abs(A[i][k]) > max_val:
                max_val = abs(A[i][k])
                max_row = i
                
        # swap rows k and max_row
        if max_row != k:
            A[k], A[max_row] = A[max_row], A[k]
            b[k], b[max_row] = b[max_row], b[k]
            
        for i in range(k + 1, n):
            if A[k][k] == 0:
                continue
            factor = A[i][k] / A[k][k]
            for j in range(k, n):
                A[i][j] -= factor * A[k][j]
            b[i] -= factor * b[k]
            
    # Back substitution
    x_sol = [0.0] * n
    for i in range(n - 1, -1, -1):
        s = sum(A[i][j] * x_sol[j] for j in range(i + 1, n))
        x_sol[i] = (b[i] - s) / A[i][i]
        
    return x_sol

def polynomial(x_val, coef):
    return sum(c * (x_val ** i) for i, c in enumerate(coef))

def polynomial_array(x, coef):
    return [polynomial(xi, coef) for xi in x]

def variance(y_true, y_approx):
    n = len(y_true)
    return sum((yt - ya) ** 2 for yt, ya in zip(y_true, y_approx)) / n

def main():
    # 1. Read input data
    x, y = read_data('data.csv')
    
    # 2. Add some headers to make it clean
    print("Лабораторна робота №3")
    print("Знаходження алгебраїчних многочленів найкращого квадратичного наближення\n")
    
    # 3. Choose optimal polynomial degree
    max_degree = 4
    variances = []
    models = []
    
    print("M \t Variance \t Coefficients")
    print("-" * 80)
    for m in range(1, max_degree + 1):
        A = form_matrix(x, m)
        b = form_vector(x, y, m)
        coef = gauss_solve(A, b)
        y_approx = polynomial_array(x, coef)
        var = variance(y, y_approx)
        variances.append(var)
        models.append(coef)
        
        coef_str = ", ".join(f"{c:10.4f}" for c in coef)
        print(f"{m} \t {var:8.4f} \t [{coef_str}]")
        
    optimal_m = variances.index(min(variances)) + 1
    optimal_coef = models[optimal_m - 1]
    
    print("-" * 80)
    print(f"Оптимальний степінь многочлена (мінімальна дисперсія): {optimal_m}")
    
    # 4. Final approximation with optimal degree
    y_approx_opt = polynomial_array(x, optimal_coef)
    
    # 5. Prediction for next 3 months
    x_future = [25, 26, 27]
    y_future = polynomial_array(x_future, optimal_coef) # we predict based on the polynomial
    
    print(f"Прогноз на 3 місяці {x_future}: {[round(yf, 2) for yf in y_future]}")
    
    # 6. Approximation error
    error = [yt - ya for yt, ya in zip(y, y_approx_opt)]
    print("\nПохибки апроксимації для оптимального степеня (m={}):".format(optimal_m))
    for xi, err in zip(x, error):
        print(f"x={xi:>2}, error={err:>8.4f}")
        
    # 7. Plots
    plt.figure(figsize=(12, 6))
    
    # Approximation and Actual Data
    plt.subplot(1, 2, 1)
    
    # More points for a smooth curve
    x_smooth = [1 + i * 0.1 for i in range(260 + 1)] # up to 27
    y_smooth = polynomial_array(x_smooth, optimal_coef)
    
    plt.plot(x, y, 'ko', label='Фактичні дані')
    plt.plot(x_smooth, y_smooth, 'b-', label=f'Апроксимація (m={optimal_m})')
    plt.plot(x_future, y_future, 'rx', markersize=8, label='Екстраполяція (Прогноз)', markeredgewidth=2)
    plt.title('Апроксимація та фактичні дані')
    plt.xlabel('Місяць')
    plt.ylabel('Температура')
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.7)
    
    # Error
    plt.subplot(1, 2, 2)
    plt.plot(x, error, 'o-', color='red', label='Похибка')
    plt.axhline(0, color='black', linestyle='--', linewidth=1)
    plt.title('Похибка апроксимації')
    plt.xlabel('Місяць')
    plt.ylabel('Похибка ($y - P(x)$)')
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend()
    
    plt.tight_layout()
    plt.savefig('result_plot.png', dpi=300)
    print("\nГрафіки збережені у файл 'result_plot.png'.")

if __name__ == "__main__":
    main()
