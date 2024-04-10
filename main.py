from matrix import Matrix
import matplotlib.pyplot as plt
import numpy as np


def Jacobi(A,b):
    r = Matrix(1,0,0,0,1)
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
        if norm < 10**-9:
            iterations = i+1
            break
    return (r,iterations,residuum)

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    N = 10
    A = Matrix(N,12,-1,-1)
    b = Matrix(1,0,0,0)
    print(Jacobi(A,b)[0].matrix)
    print(Jacobi(A,b)[1])
    iter = list(range(Jacobi(A,b)[1]))
    plt.plot(iter,Jacobi(A,b)[2])
    plt.yscale('log')
    plt.xlabel("iteracje")
    plt.ylabel("wrtość normy residuum")
    plt.title("zmiana normy residuum w kolejnych iteracjach")

    plt.show()
# Przykładowa macierz współczynników A i wektor prawych stron b
    A1= A.to_np_array()
    b1= b.to_np_array()

# Rozwiązanie układu równań Ax = b
    x = np.linalg.solve(A1, b1)
    print("Rozwiązanie x:", x)
