import numpy as np
import matplotlib.pyplot as plt

def rosenbrock(x):
    return 100.0 * (x[1] - x[0]**2)**2 + (1.0 - x[0])**2

def system_func(x):
    f1 = x[0]**2 + x[1]**2 - 4.0
    f2 = x[1] - x[0]**2
    return f1**2 + f2**2

def hooke_jeeves(f, x0, dx0, alpha=2.0, eps_x=1e-5, max_iters=1000):
    x_base = np.array(x0, dtype=float)
    dx = np.array(dx0, dtype=float)
    
    path = [np.copy(x_base)]
    
    def explore(x_start, current_dx):
        x = np.copy(x_start)
        for i in range(len(x)):
            x[i] += current_dx[i]
            if f(x) < f(x_start):
                x_start = np.copy(x)
                continue
            x[i] -= 2 * current_dx[i]
            if f(x) < f(x_start):
                x_start = np.copy(x)
                continue
            x[i] += current_dx[i]
        return x_start

    iters = 0
    while np.max(dx) > eps_x and iters < max_iters:
        iters += 1
        x_explore = explore(x_base, dx)
        
        if f(x_explore) < f(x_base):
            # Success, try pattern search
            while True:
                x_pattern = x_base + 2.0 * (x_explore - x_base) # standard pattern move
                x_explore_pattern = explore(x_pattern, dx)
                
                if f(x_explore_pattern) < f(x_explore):
                    x_base = np.copy(x_explore)
                    x_explore = np.copy(x_explore_pattern)
                    path.append(np.copy(x_explore))
                else:
                    x_base = np.copy(x_explore)
                    path.append(np.copy(x_base))
                    break
        else:
            dx = dx / alpha
            
    return x_base, path, iters

def main():
    print("Лабораторна робота №9")
    print("\n1. Тестування методу Хука-Дживса на функції Розенброка:")
    x0_rosen = [-1.2, 1.0]
    dx0_rosen = [0.5, 0.5]
    xmin_rosen, path_rosen, iters_rosen = hooke_jeeves(rosenbrock, x0_rosen, dx0_rosen)
    print(f"Початкова точка: {x0_rosen}")
    print(f"Точка мінімуму: {xmin_rosen}")
    print(f"Значення функції: {rosenbrock(xmin_rosen)}")
    print(f"Кількість ітерацій: {iters_rosen}")

    print("\n2. Розв'язок системи нелінійних рівнянь:")
    # x^2 + y^2 = 4
    # y = x^2
    x0_sys = [1.0, 1.0]
    dx0_sys = [0.5, 0.5]
    xmin_sys, path_sys, iters_sys = hooke_jeeves(system_func, x0_sys, dx0_sys)
    print(f"Початкова точка: {x0_sys}")
    print(f"Точка мінімуму (розв'язок): {xmin_sys}")
    print(f"Значення цільової функції: {system_func(xmin_sys)}")
    print(f"Кількість ітерацій: {iters_sys}")
    
    with open("path.txt", "w") as f:
        f.write("Траєкторія спуску:\n")
        for i, p in enumerate(path_sys):
            f.write(f"Крок {i}: x = {p[0]:.6f}, y = {p[1]:.6f}\n")
    print("Траєкторію спуску збережено у файл path.txt")

    # Побудова графіка
    t = np.linspace(-2.5, 2.5, 400)
    y_parabola = t**2
    
    theta = np.linspace(0, 2*np.pi, 400)
    x_circle = 2 * np.cos(theta)
    y_circle = 2 * np.sin(theta)

    plt.figure(figsize=(8, 8))
    plt.plot(t, y_parabola, label='y = x^2')
    plt.plot(x_circle, y_circle, label='x^2 + y^2 = 4')
    
    path_sys = np.array(path_sys)
    plt.plot(path_sys[:, 0], path_sys[:, 1], 'ro-', label='Траєкторія пошуку')
    plt.plot(xmin_sys[0], xmin_sys[1], 'k*', markersize=15, label='Знайдений корінь')
    
    plt.xlim(-3, 3)
    plt.ylim(-3, 3)
    plt.grid(True)
    plt.axhline(0, color='black',linewidth=0.5)
    plt.axvline(0, color='black',linewidth=0.5)
    plt.legend()
    plt.title("Розв'язок системи нелінійних рівнянь методом Хука-Дживса")
    plt.savefig("plot.png")
    print("Графік збережено у файл plot.png")

if __name__ == "__main__":
    main()
