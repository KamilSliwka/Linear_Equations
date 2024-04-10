
from vector import Vector
class Matrix:

    def __init__(self,a1,a2,a3):
        self.n = 7
        self.a1 = a1
        self.a2 = a2
        self.a3 = a3
        self.matrix = []
        self.createMatrix()

    def createMatrix(self):
        for i in range(self.n):
            row = []
            for j in range(self.n):
                if i == j:
                    row.append(self.a1)
                elif i + 1 == j or i - 1 == j:
                    row.append(self.a2)
                elif i + 2 == j or i - 2 == j:
                    row.append(self.a3)
                else:
                    row.append(0)
            self.matrix.append(row)
    #macierz razy wektor
    def multiplicationByVector(self,vector):
            newVector = Vector(1)
            for j in range(self.n):
                element = 0
                for i in range(self.n):
                    element += vector.vector[i][0]*self.matrix[j][i]
                newVector.vector[j][0] = element
            return newVector

    def multiplication(self, matrix2):
        newMatrix = Matrix(0, 0, 0)
        for i in range(self.n):
            for j in range(matrix2.n):
                element = 0
                for k in range(self.n):
                    element += self.matrix[i][k] * matrix2.matrix[k][j]
                newMatrix.matrix[i][j] = element
        return newMatrix
    def diagonal(self):
        newDiagonal = Matrix(self.a1,0,0)
        return newDiagonal
    def upper(self):
        newUpper = Matrix(0,0,0)
        for i in range(self.n):
            for j in range(self.n):
                if i>j:
                    newUpper.matrix[i][j] = self.matrix[i][j]

        return newUpper

    def lower(self):
        newLower = Matrix(0,0,0)
        for i in range(self.n):
            for j in range(self.n):
                if i<j:
                    newLower.matrix[i][j] = self.matrix[i][j]

        return newLower

    # macierz razy Liczba
    def multiplicationByNumber(self, number):
        newMatrix = Matrix(0,0,0)
        for i in range(self.n):
            for j in range(self.n):
                    newMatrix.matrix[i][j] = self.matrix[i][j]*number
        return newMatrix

    def addition(self, matrix2):
        newMatrix = Matrix(0, 0, 0)
        for i in range(self.n):
            for j in range(self.n):
                    newMatrix.matrix[i][j] = self.matrix[i][j] + matrix2.matrix[i][j]
        return newMatrix
    def invDiag(self):
        newMatrix = Matrix(0, 0, 0)
        for i in range(self.n):
            for j in range(self.n):
                if j == i:
                    newMatrix.matrix[i][j] = 1/self.matrix[i][j]
        return newMatrix
    def isDiagonalMatrix(self):
        for i in range(self.n):
            for j in range(self.n):
                if j != i and self.matrix[i][j]!=0:
                    return False
        return True