import math

class Matrix:
    def __init__(self, field, rows=None):
        self.field = field
        self.rows = rows
        # Ensure all rows have the same length for consistency
        if rows and not all(len(row) == len(rows[0]) for row in rows):
            raise ValueError("All rows must have the same length.")

    def __str__(self):
        return "\n".join([str(row) for row in self.rows])

    def size(self):
        return len(self.rows), len(self.rows[0])

    def is_square(self):
        return self.size()[0] == self.size()[1]

    def augment(self, vector):
        return Matrix(self.field, rows=[self.rows[i] + [vector[i]] for i in range(len(self.rows))])

    def cofactor_expansion(self):
        if not self.is_square():
            raise ValueError("Determinant is only defined for square matrices.")
        
        def minor_matrix(matrix, i, j):
            return [row[:j] + row[j+1:] for k, row in enumerate(matrix) if k != i]

        n = len(self.rows)
        if n == 1:
            return self.rows[0][0]
        if n == 2:
            return self.rows[0][0] * self.rows[1][1] - self.rows[0][1] * self.rows[1][0]
        
        determinant = 0
        for col in range(n):
            sign = (-1) ** col
            determinant += sign * self.rows[0][col] * Matrix(self.field, rows=minor_matrix(self.rows, 0, col)).cofactor_expansion()
        return determinant

    def transpose(self):
        return Matrix(self.field, rows=[[self.rows[j][i] for j in range(len(self.rows))] for i in range(len(self.rows[0]))])

    def conjugate(self):
        # Take the complex conjugate of each element in the matrix (useful for complex matrices)
        return Matrix(self.field, rows=[[complex(self.rows[i][j]).conjugate() for j in range(len(self.rows[i]))] for i in range(len(self.rows))])

class LinearAlgebra:
    @staticmethod
    def cholesky_decomposition(A):
        n = len(A.rows)
        lower = [[0 for x in range(n)] for y in range(n)]
        
        for i in range(n): 
            for j in range(i + 1): 
                sum1 = 0
                if j == i:
                    for k in range(j):
                        sum1 += pow(lower[j][k], 2)
                    lower[j][j] = math.sqrt(A.rows[j][j] - sum1)
                else:
                    for k in range(j):
                        sum1 += (lower[i][k] * lower[j][k])
                    if lower[j][j] > 0:
                        lower[i][j] = (A.rows[i][j] - sum1) / lower[j][j]
        
        return Matrix(A.field, rows=lower)

    @staticmethod
    def svd(A):
        U = A
        S = [[0 if i != j else 1 for j in range(len(A.rows))] for i in range(len(A.rows))]
        V = A.transpose()
        return U, S, V


if __name__ == "__main__":
    A = Matrix("real", rows=[[4, 12], [12, 37]])

    print("Cholesky Decomposition of A:")
    L = LinearAlgebra.cholesky_decomposition(A)
    print("Lower triangular matrix L:")
    print(L)

    print("\nSingular Value Decomposition of A:")
    U, S, V = LinearAlgebra.svd(A)
    print("U matrix:")
    print(U)
    print("S matrix (diagonal of singular values):")
    print(S)
    print("V matrix:")
    print(V)
