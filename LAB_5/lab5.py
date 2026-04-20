import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt

def f(x):
    return 50 + 20 * np.sin(np.pi * x / 12) + 5 * np.exp(-0.2 * (x - 12)**2)

def simpson_composite(func, a, b, n):
    if n % 2 != 0:
        n += 1  # n must be even for Simpson's Rule
    h = (b - a) / n
    x = np.linspace(a, b, n + 1)
    y = func(x)
    S = y[0] + y[-1] + 4 * np.sum(y[1:-1:2]) + 2 * np.sum(y[2:-2:2])
    return (h / 3) * S

def adaptive_simpson(func, a, b, eps, calc_count=0):
    c = (a + b) / 2
    h = (b - a) / 2
    fa, fb, fc = func(a), func(b), func(c)
    S = (h / 3) * (fa + 4 * fc + fb)
    calc_count[0] += 3
    
    def recursive_simpson(a, b, eps, S, fa, fb, fc):
        c = (a + b) / 2
        h = (b - a) / 2
        d = (a + c) / 2
        e = (c + b) / 2
        
        fd, fe = func(d), func(e)
        calc_count[0] += 2
        
        S_left = (h / 6) * (fa + 4 * fd + fc)
        S_right = (h / 6) * (fc + 4 * fe + fb)
        S2 = S_left + S_right
        
        if abs(S2 - S) <= 15 * eps:
            return S2 + (S2 - S) / 15
        else:
            return recursive_simpson(a, c, eps / 2, S_left, fa, fc, fd) + \
                   recursive_simpson(c, b, eps / 2, S_right, fc, fb, fe)
                   
    return recursive_simpson(a, b, eps, S, fa, fb, fc)

def main():
    a, b = 0, 24
    
    exact_val, exact_err = integrate.quad(f, a, b, epsabs=1e-14, epsrel=1e-14)
    print(f"Точне значення інтеграла: {exact_val:.14f}")
    
    N_values = [2**i for i in range(1, 15)]
    errors = []
    
    print("\nЗалежність похибки від кількості розбиттів N:")
    print("N \t\t S(N) \t\t\t Похибка")
    print("-" * 65)
    
    for N in N_values:
        S_N = simpson_composite(f, a, b, N)
        err = abs(S_N - exact_val)
        errors.append(max(err, 1e-16))
        print(f"{N} \t\t {S_N:.12f} \t {err:.12e}")
    
    optimal_idx = np.argmin(errors)
    optimal_N = N_values[optimal_idx]
    
    print(f"\nОптимальне N (найменша похибка): N = {optimal_N}, Похибка = {errors[optimal_idx]:.12e}")
    
    plt.figure(figsize=(10, 6))
    plt.loglog(N_values, errors, marker='o', linestyle='-', color='r')
    plt.title("Залежність похибки методу Сімпсона\nвід кількості розбиттів N")
    plt.xlabel("N (логарифмічна шкала)")
    plt.ylabel("Похибка |S(N) - exact| (логарифмічна шкала)")
    plt.grid(True, which="both", ls="--")
    plt.tight_layout()
    plt.savefig("error_vs_N_lab5.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    N_base = 32
    print(f"\nПриймаємо N = {N_base} для застосування методів підвищення точності.")
    S_h = simpson_composite(f, a, b, N_base)
    S_h_half = simpson_composite(f, a, b, N_base * 2)
    
    err_h = abs(S_h - exact_val)
    print(f"S({N_base}) = {S_h:.12f}, Похибка: {err_h:.12e}")
    print(f"S({N_base*2}) = {S_h_half:.12f}, Похибка: {abs(S_h_half - exact_val):.12e}")
    
    S_rr = S_h_half + (S_h_half - S_h) / 15
    err_rr = abs(S_rr - exact_val)
    print("\nМетод Рунге-Ромберга:")
    print(f"Уточнене значення: {S_rr:.14f}")
    print(f"Похибка: {err_rr:.12e}")
    if err_rr > 0:
        print(f"Похибка зменшилася у {err_h / err_rr:.2f} разів")
        
    S_h_quarter = simpson_composite(f, a, b, N_base * 4)
    den = S_h_quarter - 2 * S_h_half + S_h
    if den != 0:
        S_aitken = S_h_quarter - (S_h_quarter - S_h_half)**2 / den
        err_aitken = abs(S_aitken - exact_val)
        
        ratio = (S_h_half - S_h) / (S_h_quarter - S_h_half)
        p = np.log2(ratio) if ratio > 0 else float('nan')
        
        print("\nМетод Ейткена:")
        print(f"Уточнене значення: {S_aitken:.14f}")
        print(f"Похибка: {err_aitken:.12e}")
        print(f"Оцінка порядку формули p: {p:.2f} (очікується ~4)")
        
    eps_target = 1e-10
    calc_count = [0]
    S_adaptive = adaptive_simpson(f, a, b, eps_target, calc_count)
    err_adaptive = abs(S_adaptive - exact_val)
    
    print("\nАдаптивний метод Сімпсона:")
    print(f"Задана точність eps = {eps_target}")
    print(f"Уточнене значення: {S_adaptive:.14f}")
    print(f"Досягнута похибка: {err_adaptive:.12e}")
    print(f"Кількість обчислень функції f(x): {calc_count[0]}")

if __name__ == "__main__":
    main()
