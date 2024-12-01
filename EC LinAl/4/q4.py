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

    def rank(self):
        # getting the rank by calcu the rref.
        matrix = [row[:] for row in self.rows]
        rows, cols = len(matrix), len(matrix[0])
        rank = 0
        for col in range(cols):
            pivot_row = -1
            for row in range(rank, rows):
                if matrix[row][col] != 0:
                    pivot_row = row
                    break
            if pivot_row == -1:
                continue
            matrix[rank], matrix[pivot_row] = matrix[pivot_row], matrix[rank]
            pivot = matrix[rank][col]
            matrix[rank] = [val / pivot for val in matrix[rank]]
            for row in range(rank + 1, rows):
                factor = matrix[row][col]
                matrix[row] = [matrix[row][j] - factor * matrix[rank][j] for j in range(cols)]
            rank += 1
        return rank

    def rref(self):
        # Convert the matrix to rref.
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


class LinearSystem:
    def __init__(self, A, b):
        if len(A.rows) != len(b.rows):
            raise ValueError("Matrix A and vector b must have compatible sizes.")
        self.A = A
        self.b = b

    def is_consistent(self):
        # Check if the system is consistent by comparing the ranks of A and the augmented matrix.
        augmented_matrix = Matrix("real", rows=[row + [self.b.rows[i][0]] for i, row in enumerate(self.A.rows)])
        rank_A = self.A.rank()
        rank_augmented = augmented_matrix.rank()
        return rank_A == rank_augmented

    def is_consistent_solution(self):
        if not self.is_consistent():
            raise ValueError("The system of equations is inconsistent.")
        augmented_matrix = Matrix("real", rows=[row + [self.b.rows[i][0]] for i, row in enumerate(self.A.rows)])
        rref = augmented_matrix.rref()
        n = len(rref[0]) - 1
        solution = [row[-1] for row in rref if any(row[:-1])]
        return solution if len(solution) == n else None

    def solution_set(self):
        if not self.is_consistent():
            raise ValueError("The system of equations is inconsistent.")
        augmented_matrix = Matrix("real", rows=[row + [self.b.rows[i][0]] for i, row in enumerate(self.A.rows)])
        rref = augmented_matrix.rref()
        free_vars = []
        solution_set = []
        for i, row in enumerate(rref):
            if all(value == 0 for value in row[:-1]):
                continue
            pivot_col = next((j for j, value in enumerate(row[:-1]) if value != 0), None)
            if pivot_col is not None:
                solution_set.append((pivot_col, row[-1]))
            else:
                free_vars.append(i)
        return solution_set, free_vars

    def solnPLU(self):
        if not self.is_consistent():
            raise ValueError("The system of equations is inconsistent.")
        L, U = self.A.lu_decomposition()
        y = self._forward_substitution(L, self.b)
        x = self._backward_substitution(U, y)
        return x

    @staticmethod
    def _forward_substitution(L, b):
        n = len(L.rows)
        y = [0] * n
        for i in range(n):
            y[i] = b.rows[i][0] - sum(L.rows[i][j] * y[j] for j in range(i))
        return Matrix("real", rows=[[val] for val in y])

    @staticmethod
    def _backward_substitution(U, y):
        n = len(U.rows)
        x = [0] * n
        for i in range(n - 1, -1, -1):
            x[i] = (y.rows[i][0] - sum(U.rows[i][j] * x[j] for j in range(i + 1, n))) / U.rows[i][i]
        return Matrix("real", rows=[[val] for val in x])
