
from vector import Vector
import math
import numpy as np
class Matrix:

    def __init__(self,col,a1,a2,a3,empty=0):
        self.n = 10
        self.col =col
        self.a1 = a1
        self.a2 = a2
        self.a3 = a3
        self.matrix = []
        if col == 1:
            if empty != 0:
                # initial residual vector
                self.matrix = [1 for _ in range(self.n)]
            else:
                for i in range(1,self.n+1):
                    self.matrix.append(math.sin(i * (3 + 1)))
            col = [[element] for element in self.matrix]
            self.matrix = col
        else:
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
        newMatrix = Matrix(matrix2.col,0, 0, 0,1)
        for i in range(self.n):
            for j in range(matrix2.col):
                element = 0
                for k in range(self.col):
                    element += self.matrix[i][k] * matrix2.matrix[k][j]
                newMatrix.matrix[i][j] = element
        return newMatrix
    def diagonal(self):
        newDiagonal = Matrix(self.n,self.a1,0,0)
        return newDiagonal
    def upper(self):
        newUpper = Matrix(self.n,0,0,0)
        for i in range(self.n):
            for j in range(self.n):
                if i<j:
                    newUpper.matrix[i][j] = self.matrix[i][j]

        return newUpper

    def lower(self):
        newLower = Matrix(self.n,0,0,0)
        for i in range(self.n):
            for j in range(self.n):
                if i>j:
                    newLower.matrix[i][j] = self.matrix[i][j]

        return newLower

    # macierz razy Liczba
    def multiplicationByNumber(self, number):
        newMatrix = Matrix(self.col,0,0,0)
        for i in range(self.n):
            for j in range(self.col):
                    newMatrix.matrix[i][j] = self.matrix[i][j]*number
        return newMatrix

    def addition(self, matrix2):
        newMatrix = Matrix(self.col,0, 0, 0)
        for i in range(self.n):
            for j in range(self.col):
                    newMatrix.matrix[i][j] = self.matrix[i][j] + matrix2.matrix[i][j]
        return newMatrix
    def invDiag(self):
        newMatrix = Matrix(self.n,0, 0, 0)
        for i in range(self.n):
            for j in range(self.n):
                if j == i:
                    newMatrix.matrix[i][j] = 1/self.matrix[i][j]
        return newMatrix
    def isDiagonalMatrix(self):
        for i in range(self.n):
            for j in range(self.col):
                if j != i and self.matrix[i][j]!=0:
                    return False
        return True

    def norm(self):
        value = 0
        for row in self.matrix:
            for element in row:
                try:
                    value += element ** 2
                except OverflowError as e:

                    return  math.sqrt(value)


        return math.sqrt(value)

    def to_np_array(self):
        return np.array(self.matrix)