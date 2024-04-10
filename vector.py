import math


class Vector:
    def __init__(self,empty=0):
        self.n = 7
        self.vector = []
        if empty != 0:
            #initial residual vector
            self.vector = [1 for _ in range(self.n)]
        else:
            for i in range(self.n):
                self.vector.append(math.sin(i*(3+1)))
        col = [[element] for element in self.vector]
        self.vector = col

    def addition(self, vector2):
        newVector = Vector()
        for i in range(self.n):
            newVector.vector[i][0] = self.vector[i][0] + vector2.vector[i][0]
        return newVector
    def subtraction(self, vector2):
        newVector = Vector()
        for i in range(self.n):
            newVector.vector[i][0] = self.vector[i][0] - vector2.vector[i][0]
        return newVector
    def norm(self,norm):
        for i in range(self.n):
            if self.vector[i][0]>norm:
                return False
        return True