import numpy as np
import os

def generate_diagonally_dominant_matrix(n=100):
    A = np.random.rand(n, n) * 10
    # Зробити діагональне переважання:
    # Елемент на діагоналі має бути більшим за суму абсолютних значень інших елементів рядка
    for i in range(n):
        row_sum = np.sum(np.abs(A[i, :])) - np.abs(A[i, i])
        A[i, i] = row_sum + np.random.uniform(100.0, 200.0)
    return A

def generate_and_save_data(n=100, a_filename="matrix_A.txt", b_filename="vector_B.txt"):
    A = generate_diagonally_dominant_matrix(n)
    X_exact = np.full(n, 2.5)
    B = np.dot(A, X_exact)
    
    np.savetxt(a_filename, A, fmt='%15.8f')
    np.savetxt(b_filename, B, fmt='%15.8f')
    return A, B, X_exact

def simple_iteration_method(A, B, eps=1e-14, max_iter=1000):
    n = len(B)
    X = np.full(n, 1.0)
    
    D = np.diag(np.diag(A))
    D_inv = np.diag(1 / np.diag(A))
    
    C = np.eye(n) - np.dot(D_inv, A)
    d = np.dot(D_inv, B)
    
    iterations = 0
    while iterations < max_iter:
        X_new = np.dot(C, X) + d
        
        diff = np.max(np.abs(X_new - X))
        R = np.max(np.abs(np.dot(A, X_new) - B))
        
        if diff <= eps:
            return X_new, iterations + 1
            
        X = X_new
        iterations += 1
        
    return X, iterations

def jacobi_method(A, B, eps=1e-14, max_iter=1000):
    n = len(B)
    X = np.full(n, 1.0)
    X_new = np.zeros(n)
    
    iterations = 0
    while iterations < max_iter:
        for i in range(n):
            s = 0
            for j in range(n):
                if i != j:
                    s += A[i, j] * X[j]
            X_new[i] = (B[i] - s) / A[i, i]
            
        diff = np.max(np.abs(X_new - X))
        R = np.max(np.abs(np.dot(A, X_new) - B))
        
        if diff <= eps:
            return X_new, iterations + 1
            
        X = np.copy(X_new)
        iterations += 1
        
    return X, iterations

def gauss_seidel_method(A, B, eps=1e-14, max_iter=1000):
    n = len(B)
    X = np.full(n, 1.0)
    
    iterations = 0
    while iterations < max_iter:
        X_prev = np.copy(X)
        for i in range(n):
            s1 = sum(A[i, j] * X[j] for j in range(i))
            s2 = sum(A[i, j] * X_prev[j] for j in range(i + 1, n))
            X[i] = (B[i] - s1 - s2) / A[i, i]
            
        diff = np.max(np.abs(X - X_prev))
        R = np.max(np.abs(np.dot(A, X) - B))
        
        if diff <= eps:
            return X, iterations + 1
            
        iterations += 1
        
    return X, iterations

def main():
    n = 100
    print(f"Генерація матриці з діагональним переважанням {n}x{n}...\n")
    A, B, X_exact = generate_and_save_data(n)
    
    print("----- Метод простої ітерації -----")
    X_si, iters_si = simple_iteration_method(A, B, eps=1e-14)
    eps_si = np.max(np.abs(np.dot(A, X_si) - B))
    print(f"Кількість ітерацій: {iters_si}")
    print(f"Максимальна нев'язка: {eps_si:.4e}")
    print(f"Похибка відносно точного розв'язку: {np.max(np.abs(X_si - X_exact)):.4e}\n")
    
    print("----- Метод Якобі -----")
    X_jacobi, iters_jacobi = jacobi_method(A, B, eps=1e-14)
    eps_jacobi = np.max(np.abs(np.dot(A, X_jacobi) - B))
    print(f"Кількість ітерацій: {iters_jacobi}")
    print(f"Максимальна нев'язка: {eps_jacobi:.4e}")
    print(f"Похибка відносно точного розв'язку: {np.max(np.abs(X_jacobi - X_exact)):.4e}\n")
    
    print("----- Метод Зейделя -----")
    X_seidel, iters_seidel = gauss_seidel_method(A, B, eps=1e-14)
    eps_seidel = np.max(np.abs(np.dot(A, X_seidel) - B))
    print(f"Кількість ітерацій: {iters_seidel}")
    print(f"Максимальна нев'язка: {eps_seidel:.4e}")
    print(f"Похибка відносно точного розв'язку: {np.max(np.abs(X_seidel - X_exact)):.4e}\n")

    print("\nВисновок: ")
    print("Метод Зейделя використовує щойно оновлені значення змінних для швидшої збіжності, що має демонструвати меншу кількість ітерацій порівняно з Якобі.")

if __name__ == "__main__":
    main()
