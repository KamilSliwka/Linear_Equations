from matrix import Matrix
import matplotlib.pyplot as plt
import numpy as np
import time


def Jacobi(A, b):
    r = Matrix(1, 0, 0, 0, 1)
    U = A.upper()
    L = A.lower()
    D = A.diagonal()
    D = D.invDiag()
    M = U.addition(L).multiplication(D.multiplicationByNumber(-1))
    bm = D.multiplication(b)
    residuum = []
    iterations = 1000
    for i in range(iterations):
        r = M.multiplication(r)
        r = bm.addition(r)
        result = A.multiplication(r)
        errorNorm = result.addition(b.multiplicationByNumber(-1))
        norm = errorNorm.norm()
        residuum.append(norm)
        if norm < 10 ** -9:
            iterations = i + 1
            break
    return (r, iterations, residuum)


def GaussSeidl(A, b):
    r = Matrix(1, 0, 0, 0, 1)
    U = A.upper()
    L = A.lower()
    D = A.diagonal()
    L1 = D.addition(L)
    M = forwardSubstitution(L1, U)
    M = M.multiplicationByNumber(-1)
    bm = forwardSubstitution(L1, b)
    residuum = []
    iterations = 1000
    for i in range(iterations):
        r = M.multiplication(r)
        r = bm.addition(r)
        result = A.multiplication(r)
        errorNorm = result.addition(b.multiplicationByNumber(-1))
        norm = errorNorm.norm()
        residuum.append(norm)
        if norm < 10 ** -9:
            iterations = i + 1
            break
    return (r, iterations, residuum)

def distributionLU(A,b):
    L = Matrix(A.n, 1, 0, 0)
    U = A
    for k in range(A.n-1):
        for j in range(k+1,A.n):
            L.matrix[j][k] = U.matrix[j][k]/U.matrix[k][k]
            #U.matrix[j][k:A.n] = U.matrix[j][k:A.n] - L.matrix[j][k]*U.matrix[j][k:A.n]
            # multiplier = L.matrix[j][k]
            # sequence = U.matrix[j][k:A.n]
            # U.matrix[j][k:A.n] =[x-x * multiplier for x in sequence]
            for m in range(k,A.n):
                U.matrix[j][m] = U.matrix[j][m] - U.matrix[j][m]*L.matrix[j][k]
    return U,L
def factorizationLU(A,b):
    [U,L] = distributionLU(A,b)
    Y = forwardSubstitution(L,b)
    X = backwardSubstitution(U,Y)
    result = A.multiplication(X)
    errorNorm = result.addition(b.multiplicationByNumber(-1))
    norm = errorNorm.norm()

    return X,norm

def forwardSubstitution(L, B):
    r = Matrix(B.col, 0, 0, 0, 1)
    for k in range(B.col):
        for i in range(L.n):
            x = 0
            for j in range(i):
                x += L.matrix[i][j] * r.matrix[j][k]
            r.matrix[i][k] = (B.matrix[i][k] - x) / L.matrix[i][i]
    return r

def backwardSubstitution(U, B):
    r = Matrix(B.col, 0, 0, 0, 1)
    for k in range(B.col):
        for i in range(U.n-1,-1,-1):
            x = 0
            for j in range(U.n-1,i,-1):
                x += U.matrix[i][j] * r.matrix[j][k]
            r.matrix[i][k] = (B.matrix[i][k] - x) / U.matrix[i][i]
    return r




# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    N = 10
    A = Matrix(N, 3, -1, -1)
    b = Matrix(1, 0, 0, 0)
    start_time = time.time()
    [r, iterations, residuum] = Jacobi(A, b)
    end_time = time.time()
    execution_time = end_time - start_time
    print("Czas obliczeń metodą Jacobiego: ",execution_time," ")

    print(r.matrix)
    print(iterations)
    iter = list(range(iterations))
    plt.plot(iter, residuum)
    plt.yscale('log')
    plt.xlabel("iteracje")
    plt.ylabel("wrtość normy residuum")
    plt.title("zmiana normy residuum w kolejnych iteracjach [Jacobi]")
    plt.xticks(range(0, len(iter) + 1, 1))
    plt.show()

    start_time = time.time()
    [r, iterations, residuum] = GaussSeidl(A, b)
    end_time = time.time()
    execution_time = end_time - start_time
    print("Czas obliczeń metodą Gaussa-Seidla: ", execution_time, " ")

    print(r.matrix)
    print(iterations)
    iter = list(range(iterations))
    plt.plot(iter, residuum)
    plt.yscale('log')
    plt.xlabel("iteracje")
    plt.ylabel("wrtość normy residuum")
    plt.title("zmiana normy residuum w kolejnych iteracjach [Gauss-Seidl]")
    plt.xticks(range(0, len(iter) + 1, 1))
    plt.show()
    # Przykładowa macierz współczynników A i wektor prawych stron b

    A1 = A.to_np_array()
    b1 = b.to_np_array()

    # Rozwiązanie układu równań Ax = b
    x = np.linalg.solve(A1, b1)
    print("Rozwiązanie x:", x)

    start_time = time.time()
    [r, residuum] = factorizationLU(A, b)
    end_time = time.time()
    execution_time = end_time - start_time
    print("LU: ",r.matrix)
    print(residuum)
    print("Czas obliczeń metodą faktoryzacji LU: ", execution_time, " ")
