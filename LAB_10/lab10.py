import numpy as np
import matplotlib.pyplot as plt

def f(x, y):
    return -y + x + 1

def y_exact(x):
    return np.exp(-x) + x

def rk4_step(x, y, h):
    k1 = f(x, y)
    k2 = f(x + h/2, y + h/2 * k1)
    k3 = f(x + h/2, y + h/2 * k2)
    k4 = f(x + h, y + h * k3)
    return y + h/6 * (k1 + 2*k2 + 2*k3 + k4)

# ================= PART 1: Adams Predictor-Corrector =================
def adams_pc_fixed_step(a, b, y0, h):
    xs = np.arange(a, b + h, h)
    n = len(xs)
    ys = np.zeros(n)
    ys[0] = y0
    
    if n > 1:
        # Use exact or RK4 for the first step
        ys[1] = rk4_step(xs[0], ys[0], h)
        
    y_pred = np.zeros(n)
    y_pred[0], y_pred[1] = ys[0], ys[1]
    
    for i in range(1, n - 1):
        # Predictor
        fn = f(xs[i], ys[i])
        fn_minus_1 = f(xs[i-1], ys[i-1])
        yp = ys[i] + h/2 * (3*fn - fn_minus_1)
        y_pred[i+1] = yp
        
        # Corrector
        fn_plus_1 = f(xs[i+1], yp)
        yc = ys[i] + h/2 * (fn_plus_1 + fn)
        ys[i+1] = yc
        
    return xs, ys, y_pred

def adams_pc_auto_step(a, b, y0, h0, eps):
    xs = [a]
    ys = [y0]
    hs = [h0]
    
    h = h0
    x = a
    
    need_restart = True
    
    while x < b:
        if x + h > b:
            h = b - x
            need_restart = True
            
        if need_restart:
            # When changing step size, we use RK4 to take the first step 
            # and establish a reliable backward point for the next Adams iteration
            y_next = rk4_step(xs[-1], ys[-1], h)
            xs.append(x + h)
            ys.append(y_next)
            hs.append(h)
            x += h
            need_restart = False
            if x >= b: break
            continue
            
        # Normal Adams Predictor-Corrector step
        fn = f(xs[-1], ys[-1])
        fn_minus_1 = f(xs[-2], ys[-2])
        
        yp = ys[-1] + h/2 * (3*fn - fn_minus_1)
        fn_plus_1 = f(x + h, yp)
        yc = ys[-1] + h/2 * (fn_plus_1 + fn)
        
        # Error estimate
        err = abs(yc - yp) / 6.0
        
        if err > eps:
            h /= 2
            need_restart = True
            continue # Retry from current x with halved step
            
        # Accept step
        xs.append(x + h)
        ys.append(yc)
        hs.append(h)
        x += h
        
        if err < eps / 4: 
            h *= 2
            need_restart = True
            
    return np.array(xs), np.array(ys), np.array(hs)

# ================= PART 2: Runge-Kutta 4 =================
def rk4_fixed_step(a, b, y0, h):
    xs = np.arange(a, b + h, h)
    ys = np.zeros(len(xs))
    ys[0] = y0
    for i in range(len(xs) - 1):
        ys[i+1] = rk4_step(xs[i], ys[i], h)
    return xs, ys

def rk4_auto_step(a, b, y0, h0, eps):
    xs = [a]
    ys = [y0]
    hs = [h0]
    
    h = h0
    x = a
    y = y0
    
    while x < b:
        if x + h > b:
            h = b - x
            
        # 1 step with h
        y_h = rk4_step(x, y, h)
        
        # 2 steps with h/2
        y_h2_half = rk4_step(x, y, h/2)
        y_h2 = rk4_step(x + h/2, y_h2_half, h/2)
        
        err = abs(y_h - y_h2) / 15.0
        
        if err > eps:
            h /= 2
            continue
            
        x += h
        y = y_h2 # more accurate
        xs.append(x)
        ys.append(y)
        hs.append(h)
        
        if err < eps / 32: # RK4 threshold to increase
            h *= 2
            
    return np.array(xs), np.array(ys), np.array(hs)

def main():
    print("Лабораторна робота №10")
    a, b = 0, 5
    y0 = 1.0
    h_fixed = 0.1
    eps = 1e-4

    # ========== PART 1 ==========
    print("\nЧ.1. Метод Адамса (прогноз-корекція 2-го порядку)")
    xs_adams, ys_adams, yp_adams = adams_pc_fixed_step(a, b, y0, h_fixed)
    y_ex_adams = y_exact(xs_adams)
    
    err_exact_adams = np.abs(ys_adams - y_ex_adams)
    err_est_adams = np.abs(ys_adams - yp_adams) / 6.0
    
    plt.figure(figsize=(10, 5))
    plt.subplot(1, 2, 1)
    plt.plot(xs_adams, err_exact_adams, label='Точна похибка')
    plt.plot(xs_adams, err_est_adams, '--', label='Оцінка (1/6|yc - yp|)')
    plt.title("Локальна похибка методу Адамса")
    plt.xlabel('x')
    plt.ylabel('Похибка')
    plt.legend()
    plt.grid(True)
    
    xs_auto_ad, ys_auto_ad, hs_auto_ad = adams_pc_auto_step(a, b, y0, h_fixed, eps)
    plt.subplot(1, 2, 2)
    plt.plot(xs_auto_ad, hs_auto_ad, 'g.-')
    plt.title("Зміна кроку h(x) (Адамс)")
    plt.xlabel('x')
    plt.ylabel('h')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("lab10_part1_adams.png")

    # ========== PART 2 ==========
    print("Ч.2. Метод Рунге-Кутта 4-го порядку")
    h_rk = 0.01
    xs_rk, ys_rk = rk4_fixed_step(a, b, y0, h_rk)
    y_ex_rk = y_exact(xs_rk)
    err_exact_rk = np.abs(ys_rk - y_ex_rk)
    
    # Runge estimate: compute ys_rk with h/2
    xs_rk_half, ys_rk_half_steps = rk4_fixed_step(a, b, y0, h_rk/2)
    ys_rk_half = ys_rk_half_steps[::2] # take every second point to match xs_rk
    err_est_rk = np.abs(ys_rk - ys_rk_half) / 15.0
    
    plt.figure(figsize=(10, 5))
    plt.subplot(1, 2, 1)
    plt.plot(xs_rk, err_exact_rk, label='Точна похибка')
    plt.plot(xs_rk, err_est_rk, '--', label='Оцінка за Рунге')
    plt.title("Локальна похибка РК4 (h=0.01)")
    plt.xlabel('x')
    plt.ylabel('Похибка')
    plt.legend()
    plt.grid(True)
    
    xs_auto_rk, ys_auto_rk, hs_auto_rk = rk4_auto_step(a, b, y0, 0.1, eps)
    plt.subplot(1, 2, 2)
    plt.plot(xs_auto_rk, hs_auto_rk, 'm.-')
    plt.title("Зміна кроку h(x) (РК4)")
    plt.xlabel('x')
    plt.ylabel('h')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("lab10_part2_rk4.png")

    print("Графіки збережено у файли lab10_part1_adams.png та lab10_part2_rk4.png")

if __name__ == "__main__":
    main()
