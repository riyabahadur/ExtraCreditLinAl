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
        return Matrix(self.field, rows=[self.rows[i] + [vector[i][0]] for i in range(len(self.rows))])

    def solve(self, vector):
        augmented = self.augment(vector)
        rref = augmented.rref()
        # Extract the last column from the reduced row-echelon form to get the solution.
        return [row[-1] for row in rref]


class DeterminantMethods:
    @staticmethod
    def det_cofactor(A):
        """Compute the determinant using cofactor expansion."""
        if not A.is_square():
            raise ValueError("Determinant is only defined for square matrices.")

        def minor_matrix(matrix, i, j):
            # Removes row i and column j to create the minor matrix.
            return [row[:j] + row[j+1:] for k, row in enumerate(matrix) if k != i]

        n = len(A.rows)
        if n == 1:
            return A.rows[0][0]
        if n == 2:
            return A.rows[0][0] * A.rows[1][1] - A.rows[0][1] * A.rows[1][0]

        determinant = 0
        for col in range(n):
            sign = (-1) ** col
            determinant += sign * A.rows[0][col] * Matrix(A.field, rows=minor_matrix(A.rows, 0, col)).det_cofactor()
        return determinant

    @staticmethod
    def det_PLU(A):
        """Compute the determinant using PLU decomposition."""
        if not A.is_square():
            raise ValueError("Determinant is only defined for square matrices.")

        L, U = A.lu_decomposition()
        det_L = 1
        det_U = 1
        for i in range(len(L.rows)):
            det_L *= L.rows[i][i]
            det_U *= U.rows[i][i]
        # The determinant is the product of the determinants of L and U.
        return det_L * det_U

    @staticmethod
    def det_RREF(A):
        """Compute the determinant using RREF and elementary matrices."""
        if not A.is_square():
            raise ValueError("Determinant is only defined for square matrices.")

        augmented = Matrix("real", rows=[row[:] for row in A.rows])
        rref_matrix = augmented.rref()

        det = 1
        row_swap_count = 0

        for i in range(len(rref_matrix)):
            if rref_matrix[i][i] == 0:
                return 0
            det *= rref_matrix[i][i]
            if i != row_swap_count:
                row_swap_count += 1
                det *= -1  # Adjust the determinant sign for each row swap
        return det

    def lu_decomposition(self):
        rows, cols = self.size()
        if rows != cols:
            raise ValueError("LU decomposition is only defined for square matrices.")

        n = rows
        L = [[0.0] * n for _ in range(n)]
        U = [[0.0] * n for _ in range(n)]

        for i in range(n):
            for j in range(i, n):
                # Calculate U[i][j] by subtracting the dot product of the previous rows.
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
            # Swap rows to ensure the pivot element is non-zero.
            matrix[i], matrix[r] = matrix[r], matrix[i]
            lv = matrix[r][lead]
            matrix[r] = [m / float(lv) for m in matrix[r]]
            for i in range(rows):
                if i != r:
                    lv = matrix[i][lead]
                    matrix[i] = [iv - lv * rv for rv, iv in zip(matrix[r], matrix[i])]
            lead += 1
        return matrix
