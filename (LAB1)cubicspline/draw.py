import matplotlib.pyplot as plt

def plot_spline_results():
    # --- 1 Зчитування даних
    
    x_input = []
    y_input = []

    try:
        with open('input.txt', 'r') as f_in:
            for line in f_in:
                parts = line.split()
                if len(parts) >= 3:
                    x_input.append(float(parts[1]))
                    y_input.append(float(parts[2]))
    except FileNotFoundError:
        print("Помилка: Файл 'input.txt' не знайдено.")
        return

    x_out = []
    y_spline = []
    y_real = []
    errors = []

    try:
        with open('output.txt', 'r') as f_out:
            for line in f_out:
                parts = line.split()
                if len(parts) >= 5:
                    x_out.append(float(parts[1]))
                    y_real.append(float(parts[2]))
                    y_spline.append(float(parts[3]))
                    errors.append(float(parts[4]))
    except FileNotFoundError:
        print("Помилка: Файл 'output.txt' не знайдено.")
        return

    # --- 2. Побудова графіків ---
    
    # Встановлюємо 'Agg' бекенд, щоб уникнути попереджень про відсутність вікна
    # Це треба зробити ДО створення фігури, якщо matplotlib ще не ініціалізовано
    # Але зазвичай import matplotlib.pyplot as plt достатньо, просто не викликаємо show()
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Графік 1: Інтерполяція
    ax1.set_title("Кубічний сплайн та вхідні точки")
    ax1.set_xlabel("x")
    ax1.set_ylabel("y")
    
    ax1.plot(x_out, y_spline, 'r-', label='Сплайн S(x)')
    ax1.plot(x_out, y_real, 'g--', label='Справжній sin(x)', alpha=0.5)
    ax1.plot(x_input, y_input, 'bo', label='Вхідні точки', markersize=4)
    
    ax1.legend()
    ax1.grid(True)

    # Графік 2: Похибка
    ax2.set_title("Похибка інтерполяції (eps)")
    ax2.set_xlabel("x")
    ax2.set_ylabel("|S(x) - sin(x)|")
    
    ax2.plot(x_out, errors, 'm-', label='Похибка')
    ax2.legend()
    ax2.grid(True)

    # --- 3. Збереження результату ---
    
    # Замість plt.show(), зберігаємо у файл
    output_filename = "spline_graph.png"
    plt.savefig(output_filename, dpi=150) # dpi=150 робить картинку якіснішою
    
    print(f"Готово! Графік успішно збережено у файл: {output_filename}")
    
    # Закриваємо фігуру, щоб звільнити пам'ять
    plt.close()

if __name__ == "__main__":
    plot_spline_results()