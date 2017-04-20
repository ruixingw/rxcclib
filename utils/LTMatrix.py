import numpy as np


class LTMatrix(list):
    def __init__(self,L):
        """
        Accept a list of elements in a lower triangular matrix.
        """
        list.__init__(self,L)
        self.list = L
        i,j = LTMatrix.getRC(len(L) - 1)
        assert i == j, "Not a LTMatrix"
        self.dimension = i + 1

    def __getitem__(self,key):
        """
        Accept one or two integers.
        ONE: get item at the given position (count from zero)
        TWO: get item at the given (row, column) (count from zero)
        """
        if type(key) is tuple:
            return self[LTMatrix.getNum(*key)]
        return self.list[key]


    @staticmethod
    def getRC(N):
        """
        Return the row and column number of the Nth entry  of a lower triangular matrix.
        N, ROW, COLUMN are counted from ZERO!
        Example:
        0
        1  2
        3  4  5
        6  7  8  9
        10 11 12 13 14
        15 16 17 18 19 20

        >>> getRC(18)
        (5, 3)

        18th element is at row 5 and column 3. (count from zero)
        """
        N += 1
        y = int((np.sqrt(1 + 8 * N) - 1) / 2)
        b = int(N - (y**2 + y) / 2)
        if b == 0:
            return (y - 1, y - 1)
        else:
            return (y, b - 1)

    @staticmethod
    def getNum(i, j):
        """
        Return the number of entry in the i-th row and j-th column of a symmetric matrix. (count from zero)

        >>> getNum(3, 4)
        13
        """
        i += 1
        j += 1
        if i < j:
            i, j = j, i
        num = (i * (i - 1) / 2) + j
        num = int(num)
        return num - 1


    def constfullmat(self):
        L = []
        for i in range(self.dimension):
            L.append([])
            for j in range(self.dimension):
                L[-1].append(self[i,j])
        return np.array(L)

    @property
    def fullmat(self, renew=False):
        """
        Return the full matrix (in type of np.ndarray)
        """
        return getattr(self, '_fullmat', self.constfullmat())


    def inverse(self):
        return np.linalg.inv(self.fullmat)

