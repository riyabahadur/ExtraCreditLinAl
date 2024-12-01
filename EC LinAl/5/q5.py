class Matrix:
    def __init__(self, field, rows=None):
        # Initialize the matrix with a field type and row data.
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

    def determinant(self):
        if not self.is_square():
            raise ValueError("Determinant is only defined for square matrices.")
        n = len(self.rows)
        if n == 1:
            return self.rows[0][0]
        if n == 2:
            return self.rows[0][0] * self.rows[1][1] - self.rows[0][1] * self.rows[1][0]
        det = 0
        for col in range(n):
            # Minor matrix is formed by removing the first row and the current column.
            sub_matrix = Matrix(self.field, rows=[row[:col] + row[col+1:] for row in self.rows[1:]])
            det += ((-1) ** col) * self.rows[0][col] * sub_matrix.determinant()
        return det

    def inv(self):
        if not self.is_square():
            raise ValueError("Inverse is only defined for square matrices.")
        n = len(self.rows)
        identity = [[1 if i == j else 0 for j in range(n)] for i in range(n)]
        augmented = [self.rows[i] + identity[i] for i in range(n)]
        for i in range(n):
            if augmented[i][i] == 0:
                # Swap rows to avoid division by zero if the pivot is zero.
                for j in range(i + 1, n):
                    if augmented[j][i] != 0:
                        augmented[i], augmented[j] = augmented[j], augmented[i]
                        break
                else:
                    raise ValueError("Matrix is not invertible.")
            pivot = augmented[i][i]
            augmented[i] = [x / pivot for x in augmented[i]]
            for j in range(n):
                if i != j:
                    factor = augmented[j][i]
                    augmented[j] = [augmented[j][k] - factor * augmented[i][k] for k in range(2 * n)]
        inverse = [row[n:] for row in augmented]
        return Matrix(self.field, rows=inverse)

    def adjoint(self):
        if not self.is_square():
            raise ValueError("Adjoint is only defined for square matrices.")
        n = len(self.rows)
        adj = [[0] * n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                # Cofactor for each element is calculated using a minor matrix.
                sub_matrix = Matrix(self.field, rows=[row[:j] + row[j+1:] for row in self.rows[:i] + self.rows[i+1:]])
                adj[j][i] = ((-1) ** (i + j)) * sub_matrix.determinant()
        return Matrix(self.field, rows=adj)

    def invadj(self):
        if not self.is_square():
            raise ValueError("Inverse is only defined for square matrices.")
        det = self.determinant()
        if det == 0:
            raise ValueError("Matrix is not invertible.")
        adj = self.adjoint()
        # The inverse is calculated as adjoint divided by the determinant.
        inverse = [[adj.rows[i][j] / det for j in range(len(adj.rows[0]))] for i in range(len(adj.rows))]
        return Matrix(self.field, rows=inverse)
