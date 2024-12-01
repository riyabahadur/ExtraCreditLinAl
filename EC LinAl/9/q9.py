import cmath

class Matrix:
    def __init__(self, field, rows=None):
        self.field = field
        self.rows = rows
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

    def lu_decomposition(self):
        rows, cols = self.size()
        if rows != cols:
            raise ValueError("LU decomposition is only defined for square matrices.")
        
        n = rows
        L = [[0.0] * n for _ in range(n)]
        U = [[0.0] * n for _ in range(n)]

        for i in range(n):
            for j in range(i, n):
                U[i][j] = self.rows[i][j] - sum(L[i][k] * U[k][j] for k in range(i))
            for j in range(i, n):
                if i == j:
                    L[i][i] = 1.0
                else:
                    L[j][i] = (self.rows[j][i] - sum(L[j][k] * U[k][i] for k in range(i))) / U[i][i]

        return Matrix(self.field, rows=L), Matrix(self.field, rows=U)

    def rref(self):
        matrix = [row[:] for row in self.rows]
        rows, cols = len(matrix), len(matrix[0])
        lead = 0
        for r in range(rows):
            if lead >= cols:
                break
            i = r
            while matrix[i][lead] == 0:
                i += 1
                if i == rows:
                    i = r
                    lead += 1
                    if lead == cols:
                        break
            matrix[i], matrix[r] = matrix[r], matrix[i]
            lv = matrix[r][lead]
            matrix[r] = [m / float(lv) for m in matrix[r]]
            for i in range(rows):
                if i != r:
                    lv = matrix[i][lead]
                    matrix[i] = [iv - lv * rv for rv, iv in zip(matrix[r], matrix[i])]
            lead += 1
        return matrix


class LinearAlgebra:
    @staticmethod
    def poly_roots(p):
        n = len(p) - 1
        roots = [complex(cmath.exp(2 * cmath.pi * 1j * k / n)) for k in range(n)]
        for k in range(20):
            for i in range(n):
                f = sum(c * (roots[i] ** j) for j, c in enumerate(p))
                f_prime = sum(j * c * (roots[i] ** (j - 1)) for j, c in enumerate(p) if j > 0)
                roots[i] = roots[i] - f / f_prime
        return roots

    @staticmethod
    def char_poly(A):
        n = len(A.rows)
        identity_matrix = Matrix(A.field, rows=[[1 if i == j else 0 for j in range(n)] for i in range(n)])
        return [(A - identity_matrix * λ).cofactor_expansion() for λ in range(n)]

    @staticmethod
    def min_poly(A):
        pass

    @staticmethod
    def eigenvalues(A):
        return [0]

    @staticmethod
    def is_similar(A, B):
        pass

    @staticmethod
    def COB_similar(A, B):
        pass

    @staticmethod
    def alg_mul(A, λ):
        pass

    @staticmethod
    def geo_mul(A, λ):
        pass

    @staticmethod
    def eigen_basis(A, λ):
        pass

    @staticmethod
    def COB_diag(A):
        pass


if __name__ == "__main__":
    coefficients = [1, 0, -1]
    roots = LinearAlgebra.poly_roots(coefficients)
    print("Roots of the polynomial:", roots)

    A = Matrix("real", rows=[[1, 2], [3, 4]])
    print("Characteristic Polynomial of A:", LinearAlgebra.char_poly(A))

    B = Matrix("real", rows=[[2, 3], [4, 5]])
    print("Are A and B similar?", LinearAlgebra.is_similar(A, B))

    λ = 1
    print("Algebraic multiplicity of λ:", LinearAlgebra.alg_mul(A, λ))
    print("Geometric multiplicity of λ:", LinearAlgebra.geo_mul(A, λ))
    print("Eigen-basis of λ:", LinearAlgebra.eigen_basis(A, λ))

    print("Change of basis to diagonalize A:", LinearAlgebra.COB_diag(A))
