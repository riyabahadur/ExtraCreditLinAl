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
            # Swapping the current row with the pivot row to simplify the elimination process.
            matrix[rank], matrix[pivot_row] = matrix[pivot_row], matrix[rank]
            pivot = matrix[rank][col]
            matrix[rank] = [val / pivot for val in matrix[rank]]
            for row in range(rank + 1, rows):
                factor = matrix[row][col]
                matrix[row] = [matrix[row][j] - factor * matrix[rank][j] for j in range(cols)]
            rank += 1
        return rank

    def is_square(self):
        return self.size()[0] == self.size()[1]

    def augment(self, vector):
        return Matrix(self.field, rows=[self.rows[i] + [vector[i][0]] for i in range(len(self.rows))])

    def solve(self, vector):
        augmented = self.augment(vector)
        rref = augmented.rref()
        # Extracts the solution by reading the last column of the reduced matrix.
        return [row[-1] for row in rref]


class LinearAlgebra:
    @staticmethod
    def is_in_linear_span(S, v):
        matrix_S = Matrix("real", rows=S)
        vector_v = Matrix("real", rows=[[val] for val in v])
        augmented = matrix_S.augment(vector_v)
        # Verifies linear span by checking if rank remains unchanged after augmentation.
        return augmented.rank() == matrix_S.rank()

    @staticmethod
    def expressin(S, v):
        if not LinearAlgebra.is_in_linear_span(S, v):
            raise ValueError("Vector v is not in the span of S.")
        matrix_S = Matrix("real", rows=S)
        vector_v = Matrix("real", rows=[[val] for val in v])
        return matrix_S.solve(vector_v)

    @staticmethod
    def is_span_equal(S1, S2):
        matrix_S1 = Matrix("real", rows=S1)
        matrix_S2 = Matrix("real", rows=S2)
        return matrix_S1.rank() == matrix_S2.rank() and matrix_S1.rank() == matrix_S1.augment(matrix_S2.rows).rank()

    @staticmethod
    def coord(B, v):
        matrix_B = Matrix("real", rows=B)
        vector_v = Matrix("real", rows=[[val] for val in v])
        return matrix_B.solve(vector_v)

    @staticmethod
    def vector_from_coords(coords, B):
        # Constructs a vector from its coordinates in a given basis by computing the linear combination.
        return [sum(c * B[i][j] for i, c in enumerate(coords)) for j in range(len(B[0]))]

    @staticmethod
    def COB(B1, B2):
        matrix_B1 = Matrix("real", rows=B1)
        matrix_B2 = Matrix("real", rows=B2)
        return matrix_B1.solve(matrix_B2)

    @staticmethod
    def change_basis(v, B1, B2):
        coords_in_B1 = LinearAlgebra.coord(B1, v)
        return LinearAlgebra.vector_from_coords(coords_in_B1, B2)