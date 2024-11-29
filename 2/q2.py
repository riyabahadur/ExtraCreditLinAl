class Matrix:
    def __init__(self, field, rows=None):
        self.field = field
        self.rows = rows
        if rows and not all(len(row) == len(rows[0]) for row in rows):
            raise ValueError("All rows must have the same length.")

    def __mul__(self, other):
        if len(self.rows[0]) != len(other.rows):
            raise ValueError("Matrix dimensions do not align for multiplication.")
        result = [[sum(self.rows[i][k] * other.rows[k][j] for k in range(len(other.rows)))
                   for j in range(len(other.rows[0]))] for i in range(len(self.rows))]
        return Matrix(self.field, rows=result)

    def transpose(self):
        return Matrix(self.field, rows=[[self.rows[j][i] for j in range(len(self.rows))] for i in range(len(self.rows[0]))])

    def conjugate(self):
        if self.field != "complex":
            raise ValueError("Conjugate is only applicable for complex matrices.")
        return Matrix(self.field, rows=[[v.conjugate() if isinstance(v, ComplexNumber) else v for v in row] for row in self.rows])

    def conjugate_transpose(self):
        return self.conjugate().transpose()

    def __str__(self):
        return "\n".join([str(row) for row in self.rows])


class MatrixProperties:
    @staticmethod
    def is_zero(matrix):
        return all(all(value == 0 for value in row) for row in matrix.rows)

    @staticmethod
    def is_symmetric(matrix):
        if not MatrixProperties.is_square(matrix):
            return False
        return matrix.rows == matrix.transpose().rows

    @staticmethod
    def is_hermitian(matrix):
        if matrix.field != "complex":
            raise ValueError("Hermitian property is only applicable to complex matrices.")
        return matrix.rows == matrix.conjugate_transpose().rows

    @staticmethod
    def is_square(matrix):
        return len(matrix.rows) == len(matrix.rows[0])

    @staticmethod
    def is_orthogonal(matrix):
        if not MatrixProperties.is_square(matrix):
            return False
        identity = Matrix("real", rows=[[1 if i == j else 0 for j in range(len(matrix.rows))] for i in range(len(matrix.rows))])
        return (matrix * matrix.transpose()).rows == identity.rows

    @staticmethod
    def is_unitary(matrix):
        if matrix.field != "complex":
            raise ValueError("Unitary property is only applicable to complex matrices.")
        if not MatrixProperties.is_square(matrix):
            return False
        identity = Matrix("complex", rows=[[1 if i == j else 0 for j in range(len(matrix.rows))] for i in range(len(matrix.rows))])
        return (matrix * matrix.conjugate_transpose()).rows == identity.rows

    @staticmethod
    def determinant(matrix):
        if not MatrixProperties.is_square(matrix):
            raise ValueError("Determinant is only defined for square matrices.")
        if len(matrix.rows) == 1:
            return matrix.rows[0][0]
        det = 0
        for col in range(len(matrix.rows)):
            sub_matrix = Matrix("real", rows=[row[:col] + row[col+1:] for row in matrix.rows[1:]])
            det += ((-1) ** col) * matrix.rows[0][col] * MatrixProperties.determinant(sub_matrix)
        return det
