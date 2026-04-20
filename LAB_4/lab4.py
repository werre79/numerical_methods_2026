import numpy as np
import matplotlib.pyplot as plt
import os

RESULTS_DIR = "."

def M(t):
    return 50 * np.exp(-0.5 * t) + 5 * np.sin(t)

def dM_dt_exact(t):
    return -25 * np.exp(-0.5 * t) + 5 * np.cos(t)

def central_difference(func, t, h):
    return (func(t + h) - func(t - h)) / (2 * h)

def main():
    t0 = 1.0
    exact_val = dM_dt_exact(t0)
    print(f"Аналітичне значення похідної M'(1) = {exact_val:.10f}\n")

    h_values = [10**(-i) for i in range(1, 15)]
    errors = []
    
    print("h \t\t D(h) \t\t\t Похибка")
    print("-" * 60)
    for h in h_values:
        D_h = central_difference(M, t0, h)
        error = abs(D_h - exact_val)
        errors.append(error)
        print(f"{h:.0e} \t {D_h:.12f} \t {error:.12e}")
    
    optimal_index = np.argmin(errors)
    optimal_h = h_values[optimal_index]
    min_error = errors[optimal_index]
    print(f"\nОптимальний крок h = {optimal_h:.0e} з мінімальною похибкою = {min_error:.12e}")
    
    plt.figure(figsize=(10, 6))
    plt.loglog(h_values, errors, marker='o', linestyle='-', color='b')
    plt.title("Залежність похибки чисельного диференціювання\nвід кроку h")
    plt.xlabel("Крок h")
    plt.ylabel("Похибка |D(h) - exact|")
    plt.grid(True, which="both", ls="--")
    plt.gca().invert_xaxis()
    plt.tight_layout()
    plt.savefig("error_vs_h_lab4.png", dpi=300, bbox_inches='tight')
    plt.close()
    h_fixed = 0.01
    D_h = central_difference(M, t0, h_fixed)
    D_h_half = central_difference(M, t0, h_fixed / 2)
    
    print(f"\nПриймаємо крок h = {h_fixed}")
    print(f"D(h) = {D_h:.10f}, D(h/2) = {D_h_half:.10f}")
    
    error_fixed_h = abs(D_h - exact_val)
    print(f"Похибка при кроці h: {error_fixed_h:.10e}")
    
    D_rr = D_h_half + (D_h_half - D_h) / 3
    error_rr = abs(D_rr - exact_val)
    print(f"\nМетод Рунге-Ромберга:")
    print(f"Уточнене значення D_RR: {D_rr:.10f}")
    print(f"Похибка методу Р-Р: {error_rr:.10e}")
    if error_fixed_h > 0:
        print(f"Похибка зменшилась у {error_fixed_h / error_rr:.2f} разів порівняно з h={h_fixed}")
    else:
        print("Похибка стала нульовою (або близькою до машинного нуля).")
        
    D_h_quarter = central_difference(M, t0, h_fixed / 4)
    D_num = D_h_half * D_h_half - D_h * D_h_quarter
    D_den = 2 * D_h_half - D_h - D_h_quarter
    
    if D_den != 0:
        D_aitken = D_h_quarter - (D_h_quarter - D_h_half)**2 / (D_h_quarter - 2*D_h_half + D_h)
        error_aitken = abs(D_aitken - exact_val)
        
        ratio = (D_h_half - D_h) / (D_h_quarter - D_h_half)
        if ratio > 0:
            p = np.log2(ratio)
        else:
            p = float('nan')
            
        print(f"\nМетод Ейткена:")
        print(f"D(h/4) = {D_h_quarter:.10f}")
        print(f"Уточнене значення D_Aitken: {D_aitken:.10f}")
        print(f"Похибка методу Ейткена: {error_aitken:.10e}")
        print(f"Оцінка порядку точності p: {p:.2f}")

if __name__ == "__main__":
    main()
