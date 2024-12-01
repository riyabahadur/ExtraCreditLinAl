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
        # Checks if the matrix is square (number of rows equals number of columns).
        return self.size()[0] == self.size()[1]

    def augment(self, vector):
        # Adds a column (vector) to the right of the matrix.
        return Matrix(self.field, rows=[self.rows[i] + [vector[i]] for i in range(len(self.rows))])

    def solve(self, vector):
        # Solves the system Ax = vector using RREF.
        augmented = self.augment(vector)
        rref = augmented.rref()
        return [row[-1] for row in rref]

    def cofactor_expansion(self):
        # Computes the determinant using cofactor expansion.
        if not self.is_square():
            raise ValueError("Determinant is only defined for square matrices.")

        def minor_matrix(matrix, i, j):
            # Removes the i-th row and j-th column to create the minor matrix.
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
        # Performs LU decomposition of the matrix.
        rows, cols = self.size()
        if rows != cols:
            raise ValueError("LU decomposition is only defined for square matrices.")

        n = rows
        L = [[0.0] * n for _ in range(n)]
        U = [[0.0] * n for _ in range(n)]

        for i in range(n):
            for j in range(i, n):
                # Compute upper triangular matrix U.
                U[i][j] = self.rows[i][j] - sum(L[i][k] * U[k][j] for k in range(i))
            for j in range(i, n):
                # Compute lower triangular matrix L.
                if i == j:
                    L[i][i] = 1.0
                else:
                    L[j][i] = (self.rows[j][i] - sum(L[j][k] * U[k][i] for k in range(i))) / U[i][i]

        return Matrix(self.field, rows=L), Matrix(self.field, rows=U)

    def rref(self):
        # Computes the Reduced Row-Echelon Form (RREF) of the matrix.
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
    def inner_product(v1, v2):
        # Calculates the inner product (dot product) of two vectors.
        if len(v1) != len(v2):
            raise ValueError("Vectors must have the same length.")
        return sum(v1[i] * v2[i] for i in range(len(v1)))

    @staticmethod
    def is_ortho(v1, v2):
        # Checks if two vectors are orthogonal.
        return LinearAlgebra.inner_product(v1, v2) == 0

    @staticmethod
    def gram_schmidt(S):
        # Performs Gram-Schmidt orthogonalization on a set of vectors.
        orthogonal_vectors = []
        for v in S:
            for u in orthogonal_vectors:
                v = [v[i] - LinearAlgebra.inner_product(v, u) * u[i] for i in range(len(v))]
            orthogonal_vectors.append(v)
        return orthogonal_vectors

    @staticmethod
    def qr_factorization(A):
        # Computes the QR factorization of matrix A.
        A = [row[:] for row in A]
        m, n = len(A), len(A[0])
        Q = [[0] * n for _ in range(m)]
        R = [[0] * n for _ in range(n)]

        for j in range(n):
            v = A[j]
            for i in range(j):
                R[i][j] = LinearAlgebra.inner_product(Q[i], v)
                v = [v[k] - R[i][j] * Q[i][k] for k in range(m)]
            R[j][j] = (sum(x * x for x in v))**0.5
            Q[j] = [x / R[j][j] for x in v]
        return Q, R

    @staticmethod
    def pseudo_inverse(S):
        # Computes the Moore-Penrose Pseudoinverse of matrix S.
        S_T = list(zip(*S))
        S_T_S = [[LinearAlgebra.inner_product(S_T[i], S_T[j]) for j in range(len(S_T))] for i in range(len(S_T))]
        S_T_S_inv = LinearAlgebra.pseudo_inverse(S_T_S)
        return [[LinearAlgebra.inner_product(S_T_S_inv[i], S_T[j]) for j in range(len(S[0]))] for i in range(len(S_T))]

    @staticmethod
    def least_square(A, b):
        # Solves the least-squares problem Ax = b using pseudo-inverse.
        A_T = list(zip(*A))
        A_T_A = [[LinearAlgebra.inner_product(A_T[i], A_T[j]) for j in range(len(A_T))] for i in range(len(A_T))]
        A_T_b = [LinearAlgebra.inner_product(A_T[i], b) for i in range(len(A_T))]
        return LinearAlgebra.pseudo_inverse(A_T_A)


if __name__ == "__main__":
    v1 = [1, 2, 3]
    v2 = [4, 5, 6]
    S = [[1, 2], [2, 3], [3, 4]]
    A = [[1, 2], [3, 4], [5, 6]]
    b = [7, 8, 9]

    print("Inner product of v1 and v2:", LinearAlgebra.inner_product(v1, v2))
    print("Are v1 and v2 orthogonal?", LinearAlgebra.is_ortho(v1, v2))

    orthogonal_vectors = LinearAlgebra.gram_schmidt(S)
    print("\nGram-Schmidt orthogonalization on S:")
    for vec in orthogonal_vectors:
        print(vec)

    Q, R = LinearAlgebra.qr_factorization(A)
    print("\nQR Factorization of A:")
    print("Q matrix:")
    for row in Q:
        print(row)
    print("R matrix:")
    for row in R:
        print(row)

    pseudo_inv_A = LinearAlgebra.pseudo_inverse(A)
    print("\nMoore-Penrose Pseudoinverse of A:")
    for row in pseudo_inv_A:
        print(row)

    least_square_solution = LinearAlgebra.least_square(A, b)
    print("\nLeast Square solution to Ax = b:")
    for row in least_square_solution:
        print(row)
