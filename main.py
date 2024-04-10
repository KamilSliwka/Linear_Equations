from matrix import Matrix
from vector import Vector



def Jacobi(A,b):
    r = Vector(1)
    U = A.upper()
    L = A.lower()
    D = A.diagonal()
    D = D.invDiag()
    M = U.addition(L).multiplication(D.multiplicationByNumber(-1))
    bm = D.multiplicationByVector(b)
    for i in range(1000):
        r = M.multiplicationByVector(r)
        r = bm.addition(r)
        result  = A.multiplicationByVector(r)
        errorNorm = result.subtraction(b)
        if errorNorm.norm(10**-12):
            break
    return r

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    A = Matrix(12,-1,-1)
    b = Vector()
    print(Jacobi(A,b).vector)



