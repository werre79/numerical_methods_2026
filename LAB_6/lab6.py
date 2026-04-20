import numpy as np
import random
import os

def generate_and_save_data(n=100, x_val=2.5, a_filename="matrix_A.txt", b_filename="vector_B.txt"):
    A = np.random.rand(n, n) * 10
    
    X_exact = np.full(n, x_val)
    
    B = np.dot(A, X_exact)
    
    np.savetxt(a_filename, A, fmt='%15.8f')
    np.savetxt(b_filename, B, fmt='%15.8f')
    
    return A, B, X_exact

def load_data(a_filename="matrix_A.txt", b_filename="vector_B.txt"):
    A = np.loadtxt(a_filename)
    B = np.loadtxt(b_filename)
    return A, B

def lu_decomposition(A):
    n = len(A)
    L = np.zeros((n, n))
    U = np.zeros((n, n))
    
    np.fill_diagonal(U, 1.0)
    
    for k in range(n):
        for i in range(k, n):
            sum_l = sum(L[i, j] * U[j, k] for j in range(k))
            L[i, k] = A[i, k] - sum_l
            
        for j in range(k + 1, n):
            sum_u = sum(L[k, idx] * U[idx, j] for idx in range(k))
            if L[k, k] == 0:
                raise ValueError("Ділення на нуль в LU-розкладі.")
            U[k, j] = (A[k, j] - sum_u) / L[k, k]
            
    return L, U

def solve_lu(L, U, B):
    n = len(B)
    Z = np.zeros(n)
    X = np.zeros(n)
    
    for i in range(n):
        Z[i] = (B[i] - sum(L[i, j] * Z[j] for j in range(i))) / L[i, i]
        
    for i in range(n - 1, -1, -1):
        X[i] = Z[i] - sum(U[i, j] * X[j] for j in range(i + 1, n))
        
    return X

def iterative_refinement(A, L, U, B, X0, eps=1e-14, max_iter=200):
    k = 0
    X = np.copy(X0)
    
    while True:
        R = B - np.dot(A, X)
        
        if np.max(np.abs(R)) <= eps:
            break
            
        dX = solve_lu(L, U, R)
        
        X_new = X + dX
        k += 1
        
        if np.max(np.abs(dX)) <= eps:
            X = X_new
            break
            
        X = X_new
        if k >= max_iter:
            print("Увага: досягнуто максимальну кількість ітерацій")
            break
            
    return X, k

def main():
    n = 100
    print(f"Генерація системи лінійних рівнянь розмірності {n}x{n}...")
    A_gen, B_gen, X_exact = generate_and_save_data(n=n, a_filename="matrix_A.txt", b_filename="vector_B.txt")
    print("Матриця A та вектор B збережені у файли matrix_A.txt, vector_B.txt")
    
    A, B = load_data("matrix_A.txt", "vector_B.txt")
    
    print("\nВиконуємо LU-розклад...")
    L, U = lu_decomposition(A)
    np.savetxt("matrix_L.txt", L, fmt='%15.8f')
    np.savetxt("matrix_U.txt", U, fmt='%15.8f')
    
    print("Розв'язуємо систему рівнянь AX = B...")
    X_lu = solve_lu(L, U, B)
    
    B_calc = np.dot(A, X_lu)
    eps_lu = np.max(np.abs(B_calc - B))
    print(f"Максимальна нев'язка (початкова) eps_lu = {eps_lu:.4e}")
    
    print("\nПроводимо ітераційне уточнення розв'язку...")
    X_refined, iterations = iterative_refinement(A, L, U, B, X_lu, eps=1e-14)
    
    B_calc_refined = np.dot(A, X_refined)
    eps_refined = np.max(np.abs(B_calc_refined - B))
    
    print(f"Кількість ітерацій: {iterations}")
    print(f"Максимальна нев'язка (після уточнення): {eps_refined:.4e}")
    
    print(f"\nПохибка відносно справжнього вектора X_exact (2.5) до уточнення: {np.max(np.abs(X_lu - X_exact)):.4e}")
    print(f"Похибка відносно справжнього вектора X_exact (2.5) після: {np.max(np.abs(X_refined - X_exact)):.4e}")

if __name__ == "__main__":
    main()
