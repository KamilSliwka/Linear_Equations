from matrix import Matrix
import matplotlib.pyplot as plt
import numpy as np
import time


def jacobi(A, b):
    r = Matrix(A.n, 1, 0, 0, 0, 1)
    iterations = 1000
    residuum = []
    for k in range(1000):
        old = r.multiplicationByNumber(1)
        for i in range(r.n):
            sum_ = 0
            for j in range(A.n):
                if j != i:
                    sum_ += A.matrix[i][j] * old.matrix[j][0]

            r.matrix[i][0] = (1 / A.matrix[i][i]) * (b.matrix[i][0] - sum_)

        result = A.multiplication(r)
        errorNorm = result.addition(b.multiplicationByNumber(-1))
        norm = errorNorm.normMaxAbs()
        residuum.append(norm)
        if norm < 10 ** -9 or norm > 10 ** 9:
            iterations = k + 1
            break
    return (r, iterations, residuum)


def Jacobi(A, b):
    r = Matrix(A.n, 1, 0, 0, 0, 1)
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
        norm = errorNorm.normMaxAbs()
        residuum.append(norm)
        if norm < 10 ** -9:
            iterations = i + 1
            break
    return (r, iterations, residuum)


def gaussSeidl(A, b):
    r = Matrix(A.n, 1, 0, 0, 0, 1)
    iterations = 1000
    residuum = []
    for k in range(1000):
        for i in range(r.n):
            sum_ = 0
            for j in range(A.n):
                if j != i:
                    sum_ += A.matrix[i][j] * r.matrix[j][0]
            r.matrix[i][0] = (1 / A.matrix[i][i]) * (b.matrix[i][0] - sum_)

        result = A.multiplication(r)
        errorNorm = result.addition(b.multiplicationByNumber(-1))
        norm = errorNorm.normMaxAbs()
        residuum.append(norm)
        if norm < 10 ** -9 or norm > 10 ** 9:
            iterations = k + 1
            break
    return (r, iterations, residuum)


def GaussSeidl(A, b):
    r = Matrix(A.n, 1, 0, 0, 0, 1)
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


def distributionLU(A):
    L = Matrix(A.n, A.n, 1, 0, 0)
    U = A.multiplicationByNumber(1)
    for k in range(A.n - 1):
        for j in range(k + 1, A.n):
            L.matrix[j][k] = U.matrix[j][k] / U.matrix[k][k]
            for m in range(k, A.n):
                U.matrix[j][m] = U.matrix[j][m] - U.matrix[k][m] * L.matrix[j][k]
    return U, L

def DistributionLU(A):
    L = Matrix(A.n, A.n, 1, 0, 0)
    U = A.multiplicationByNumber(1)
    for i in range(2,A.n+1):
        for j in range(1,i):
            L.matrix[i-1][j-1] = U.matrix[i-1][j-1]/U.matrix[j-1][j-1]
            for m in range(A.n):
                U.matrix[i-1][m] = U.matrix[i-1][m] - U.matrix[j-1][m] * L.matrix[i-1][j-1]
    return U, L


def factorizationLU(A, b):
    [U, L] = distributionLU(A)
    Y = forwardSubstitution(L, b)
    X = backwardSubstitution(U, Y)
    result = A.multiplication(X)
    errorNorm = result.addition(b.multiplicationByNumber(-1))
    norm = errorNorm.norm()

    return X, norm


def forwardSubstitution(L, B):
    r = Matrix(B.n, B.col, 0, 0, 0, 1)

    for i in range(L.n):
        x = 0
        for j in range(i):
            x += L.matrix[i][j] * r.matrix[j][0]
        r.matrix[i][0] = (B.matrix[i][0] - x) / L.matrix[i][i]

    return r


def backwardSubstitution(U, B):
    r = Matrix(B.n, B.col, 0, 0, 0, 1)
    for i in range(U.n - 1, -1, -1):
        x = 0
        for j in range(U.n - 1, i, -1):
            x += U.matrix[i][j] * r.matrix[j][0]
        r.matrix[i][0] = (B.matrix[i][0] - x) / U.matrix[i][i]
    return r


def plotTime(timeJ, timeGS, timeLU):
    plt.plot(size, timeJ, label='Metoda Jacobiego')
    plt.plot(size, timeGS, label='Metoda Gaussa-Seidla')
    plt.plot(size, timeLU, label='Metoda faktoryzacji LU')
    plt.xlabel("rozmiar macierzy")
    plt.ylabel("czas obliczeń w sekundach")
    title = "czas trwania obliczeń w zależności od rozmiaru macierzy"
    plt.title(title)
    plt.legend()
    plt.show()


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    N = 940
    A = Matrix(N, N, 3, -1, -1)
    b = Matrix(N, 1, 0, 0, 0)
    start_time = time.time()
    [r, iterationsJ, residuumJ] = jacobi(A, b)
    end_time = time.time()
    execution_time = end_time - start_time
    print("Czas obliczeń metodą Jacobiego: ", execution_time, " ")

    print(r.matrix)
    print(iterationsJ)
    iter = list(range(iterationsJ))
    plt.plot(iter, residuumJ)
    plt.yscale('log')
    plt.xlabel("iteracje")
    plt.ylabel("wrtość normy residuum")
    plt.title("zmiana normy residuum w kolejnych iteracjach [Jacobi]")
    plt.xticks(range(0, len(iter) + 1, 5))

    plt.axhline(y=1e9, color='r', linestyle='--', label='górna granica')
    plt.legend()
    plt.show()

    start_time = time.time()
    [r, iterationsG, residuumG] = gaussSeidl(A, b)
    end_time = time.time()
    execution_time = end_time - start_time
    print("Czas obliczeń metodą Gaussa-Seidla: ", execution_time, " ")

    print(r.matrix)
    print(iterationsG)
    iterG = list(range(iterationsG))
    plt.plot(iterG, residuumG)
    # plt.plot(iterG, residuumG,label='Metoda Gaussa-Seidla')
    plt.yscale('log')
    plt.xlabel("iteracje")
    plt.ylabel("wrtość normy residuum")
    plt.title("zmiana normy residuum w kolejnych iteracjach [Gauss-Seidl]")
    # plt.title("zmiana normy residuum w kolejnych iteracjach ")
    plt.xticks(range(0, len(iterG) + 1, 2))
    plt.axhline(y=1e9, color='r', linestyle='--', label='górna granica')
    plt.legend()
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
    print("LU: ", r.matrix)
    print()
    print("Wartość normy residuum: ", residuum)
    print("Czas obliczeń metodą faktoryzacji LU: ", execution_time, " s")

    size = [100, 500,1000,1500]
    timeGS = []
    timeJ = []
    timeLU = []
    for i in size:
        N = i
        A = Matrix(N, N, 12, -1, -1)
        b = Matrix(N, 1, 0, 0, 0)

        start_time = time.time()
        [r, iterations, residuum] = jacobi(A, b)
        end_time = time.time()
        execution_time = end_time - start_time
        timeJ.append(execution_time)

        start_time = time.time()
        [r, iterations, residuum] = gaussSeidl(A, b)
        end_time = time.time()
        execution_time = end_time - start_time
        timeGS.append(execution_time)

        start_time = time.time()
        [r, residuum] = factorizationLU(A, b)
        end_time = time.time()
        execution_time = end_time - start_time
        timeLU.append(execution_time)
    plotTime(timeJ, timeGS, timeLU)
