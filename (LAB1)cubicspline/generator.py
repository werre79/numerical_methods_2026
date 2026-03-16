"""
Laboratory Work #1: Cubic Spline Interpolation for Altitude Profile Analysis.
This script fetches GPS data, computes cubic splines with different node counts,
and visualizes the results along with interpolation errors.
"""
import requests
import numpy as np
import matplotlib.pyplot as plt

# Отримання даних з API
url = "https://api.open-elevation.com/api/v1/lookup?locations=48.164214,24.536044|48.164983,24.534836|48.165605,24.534068|48.166228,24.532915|48.166777,24.531927|48.167326,24.530884|48.167011,24.530061|48.166053,24.528039|48.166655,24.526064|48.166497,24.523574|48.166128,24.520214|48.165416,24.517170|48.164546,24.514640|48.163412,24.512980|48.162331,24.511715|48.162015,24.509462|48.162147,24.506932|48.161751,24.504244|48.161197,24.501793|48.160580,24.500537|48.160250,24.500106"
response = requests.get(url)
data = response.json()

results = data["results"]
n = len(results)

# Вивід вхідних координат
print("Кількість вузлів:", n)
print("\nТабуляція вузлів:")
print("№  | Latitude  | Longitude | Elevation (m)")
print("-" * 46)
for i, point in enumerate(results):
    print(f"{i:2d} | {point['latitude']:9.6f} | {point['longitude']:9.6f} | {point['elevation']:8.2f}")

# Обчислення відстаней між вузлами
def haversine(lat1, lon1, lat2, lon2):
    R = 6371000
    phi1, phi2 = np.radians(lat1), np.radians(lat2)
    dphi = np.radians(lat2 - lat1)
    dlambda = np.radians(lon2 - lon1)
    a = np.sin(dphi/2)**2 + np.cos(phi1)*np.cos(phi2)*np.sin(dlambda/2)**2
    return 2 * R * np.arctan2(np.sqrt(a), np.sqrt(1 - a))

coords = [(p["latitude"], p["longitude"]) for p in results]
elevations = [p["elevation"] for p in results]
distances = [0]
for i in range(1, n):
    d = haversine(*coords[i-1], *coords[i])
    distances.append(distances[-1] + d)

# Вивід кумулятивної відстані
print("\nТабуляція (відстань, висота):")
print("№  | Distance (m) | Elevation (m)")
print("-" * 33)
for i in range(n):
    print(f"{i:2d} | {distances[i]:12.2f} | {elevations[i]:13.2f}")

# Функція побудови кубічного сплайна
def build_spline(x_nodes, y_nodes, print_coeffs=False):
    m = len(x_nodes)
    h = [x_nodes[i] - x_nodes[i-1] for i in range(1, m)]
    
    alpha, beta, gamma, delta = np.zeros(m), np.zeros(m), np.zeros(m), np.zeros(m)
    A, B, c = np.zeros(m), np.zeros(m), np.zeros(m)
    
    # Пряма прогонка
    for i in range(1, m - 1):
        alpha[i] = h[i - 1]
        beta[i] = 2 * (h[i - 1] + h[i])
        gamma[i] = h[i]
        delta[i] = 3 * (((y_nodes[i+1] - y_nodes[i]) / h[i]) - ((y_nodes[i] - y_nodes[i-1]) / h[i-1]))
        
        A[i] = -gamma[i] / (alpha[i] * A[i-1] + beta[i])
        B[i] = (delta[i] - alpha[i] * B[i-1]) / (alpha[i] * A[i-1] + beta[i])
        
    # Зворотна прогонка
    for i in range(m - 2, 0, -1):
        c[i] = (A[i] * c[i+1]) + B[i]
        
    a_coef, b_coef, d_coef = np.zeros(m - 1), np.zeros(m - 1), np.zeros(m - 1)
    
    # Розрахунок коефіцієнтів
    for i in range(m - 1):
        a_coef[i] = y_nodes[i]
        d_coef[i] = (c[i+1] - c[i]) / (3 * h[i])
        b_coef[i] = ((y_nodes[i+1] - y_nodes[i]) / h[i]) - ((h[i] * (c[i+1] + 2 * c[i])) / 3)
        
    if print_coeffs:
        print("\nМетод прогонки (коефіцієнти матриці):")
        print(" i |   alpha  |   beta   |   gamma  |   delta  |     A    |     B    ")
        print("-" * 75)
        for i in range(1, m - 1):
            print(f"{i:2d} | {alpha[i]:8.2f} | {beta[i]:8.2f} | {gamma[i]:8.2f} | {delta[i]:8.4f} | {A[i]:8.4f} | {B[i]:8.4f}")
            
        print("\nКоефіцієнти кубічного сплайна (a, b, c, d):")
        print(" i |      a     |      b     |      c     |      d      ")
        print("-" * 65)
        for i in range(m - 1):
            print(f"{i:2d} | {a_coef[i]:10.2f} | {b_coef[i]:10.4f} | {c[i]:10.4f} | {d_coef[i]:10.6f}")

    # Генерація точок для плавного відображення
    x_smooth = np.linspace(x_nodes[0], x_nodes[-1], 1000)
    y_smooth = []
    
    for x_val in x_smooth:
        idx = np.searchsorted(x_nodes, x_val) - 1
        if idx < 0: idx = 0
        if idx >= m - 1: idx = m - 2
        
        dx = x_val - x_nodes[idx]
        y_val = a_coef[idx] + b_coef[idx]*dx + c[idx]*dx**2 + d_coef[idx]*dx**3
        y_smooth.append(y_val)
        
    return x_smooth, y_smooth

x_smooth, y_smooth = build_spline(distances, elevations, print_coeffs=True)

# Побудова графіків
# 1. Побудова окремого графіка профілю
plt.figure(figsize=(12, 7))
idx_10 = np.linspace(0, n - 1, 10, dtype=int)
x_10, y_10 = build_spline([distances[i] for i in idx_10], [elevations[i] for i in idx_10])
plt.plot(x_10, y_10, 'g--', label="Сплайн (10 вузлів)")

idx_15 = np.linspace(0, n - 1, 15, dtype=int)
x_15, y_15 = build_spline([distances[i] for i in idx_15], [elevations[i] for i in idx_15])
plt.plot(x_15, y_15, 'm-.', label="Сплайн (15 вузлів)")

idx_20 = np.linspace(0, n - 1, 20, dtype=int)
x_20, y_20 = build_spline([distances[i] for i in idx_20], [elevations[i] for i in idx_20])
plt.plot(x_20, y_20, 'b-', label="Сплайн (20 вузлів)")

plt.plot(x_smooth, y_smooth, 'k:', label="Сплайн (21 вузол - ідеальний)", linewidth=2, alpha=0.5)
plt.plot(distances, elevations, 'ro', label="Реальні GPS вузли")

plt.title("Профіль висоти: порівняння точності сплайнів")
plt.xlabel("Кумулятивна відстань (м)")
plt.ylabel("Висота (м)")
plt.legend()
plt.grid(True)
plt.savefig("goverla_profile.png", dpi=300, bbox_inches='tight')
plt.close()

# 2. Побудова окремого графіка похибок
err_10 = np.abs(np.array(y_smooth) - np.array(y_10))
err_15 = np.abs(np.array(y_smooth) - np.array(y_15))
err_20 = np.abs(np.array(y_smooth) - np.array(y_20))

plt.figure(figsize=(12, 5))
plt.plot(x_smooth, err_10, 'g--', label="Похибка (21-10 вузлів)")
plt.plot(x_smooth, err_15, 'm-.', label="Похибка (21-15 вузлів)")
plt.plot(x_smooth, err_20, 'b-', label="Похибка (21-20 вузлів)")

plt.title("Абсолютна похибка відносно 21-вузлового сплайну")
plt.xlabel("Кумулятивна відстань (м)")
plt.ylabel("Похибка (м)")
plt.legend()
plt.grid(True)
plt.savefig("goverla_errors.png", dpi=300, bbox_inches='tight')
plt.close()

# 3. Побудова комбінованого графіка (для зручності перегляду)
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10), gridspec_kw={'height_ratios': [2, 1]})

ax1.plot(x_10, y_10, 'g--', label="Сплайн (10 вузлів)")
ax1.plot(x_15, y_15, 'm-.', label="Сплайн (15 вузлів)")
ax1.plot(x_20, y_20, 'b-', label="Сплайн (20 вузлів)")
ax1.plot(x_smooth, y_smooth, 'k:', label="Сплайн (21 вузол - ідеальний)", linewidth=2, alpha=0.5)
ax1.plot(distances, elevations, 'ro', label="Реальні GPS вузли")

ax1.set_title("Профіль висоти: порівняння точності сплайнів")
ax1.set_ylabel("Висота (м)")
ax1.legend()
ax1.grid(True)

ax2.plot(x_smooth, err_10, 'g--', label="Похибка (21-10)")
ax2.plot(x_smooth, err_15, 'm-.', label="Похибка (21-15)")
ax2.plot(x_smooth, err_20, 'b-', label="Похибка (21-20)")

ax2.set_title("Абсолютна похибка відносно 21-вузлового сплайну")
ax2.set_xlabel("Кумулятивна відстань (м)")
ax2.set_ylabel("Похибка (м)")
ax2.legend()
ax2.grid(True)

plt.tight_layout()
plt.savefig("goverla_comparison.png", dpi=300, bbox_inches='tight')
print("\nГрафіки збережено (goverla_profile.png, goverla_errors.png, goverla_comparison.png)")

print("\n--- ХАРАКТЕРИСТИКИ МАРШРУТУ ---")

print(f"Загальна довжина маршруту (м): {distances[-1]:.2f}")

total_ascent = sum(max(elevations[i]-elevations[i-1], 0) for i in range(1, n))
print(f"Сумарний набір висоти (м): {total_ascent:.2f}")

total_descent = sum(max(elevations[i-1]-elevations[i], 0) for i in range(1, n))
print(f"Сумарний спуск (м): {total_descent:.2f}")

grad_full = np.gradient(y_smooth, x_smooth) * 100
print(f"Максимальний підйом (%): {np.max(grad_full):.2f}")
print(f"Максимальний спуск (%): {np.min(grad_full):.2f}")
print(f"Середній градієнт (%): {np.mean(np.abs(grad_full)):.2f}")

mass = 80
g = 9.81
energy = mass * g * total_ascent

print(f"Механічна робота (Дж): {energy:.2f}")
print(f"Механічна робота (кДж): {energy / 1000:.2f}")
print(f"Енергія (ккал): {energy / 4184:.2f}")
print("-------------------------------")