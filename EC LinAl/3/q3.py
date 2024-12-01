class ComplexNumber:
    def __init__(self, real, imag):
        self.real, self.imag = real, imag

    def __add__(self, other):
        return ComplexNumber(self.real + other.real, self.imag + other.imag)

    def __sub__(self, other):
        return ComplexNumber(self.real - other.real, self.imag - other.imag)

    def __mul__(self, other):
        return ComplexNumber(self.real * other.real - self.imag * other.imag,
                             self.real * other.imag + self.imag * other.real)

    def conjugate(self):
        return ComplexNumber(self.real, -self.imag)

    def __eq__(self, other):
        return self.real == other.real and self.imag == other.imag

    def __str__(self):
        return f"{self.real} {'+' if self.imag >= 0 else '-'} {abs(self.imag)}i"


class Matrix:
    def __init__(self, field, rows=None):
        self.field, self.rows = field, rows
        if rows and any(len(row) != len(rows[0]) for row in rows):
            raise ValueError("All rows must have the same length.")

    def __mul__(self, other):
        if len(self.rows[0]) != len(other.rows):
            raise ValueError("Matrix dimensions do not align for multiplication.")
        return Matrix(self.field, [[sum(self.rows[i][k] * other.rows[k][j] for k in range(len(other.rows)))
                                    for j in range(len(other.rows[0]))] for i in range(len(self.rows))])

    def transpose(self):
        return Matrix(self.field, [[self.rows[j][i] for j in range(len(self.rows))] for i in range(len(self.rows[0]))])

    def conjugate(self):
        if self.field != "complex":
            raise ValueError("Conjugate is only applicable for complex matrices.")
        return Matrix(self.field, [[v.conjugate() for v in row] for row in self.rows])

    def conjugate_transpose(self):
        return self.conjugate().transpose()

    def __str__(self):
        return "\n".join(map(str, self.rows))


class MatrixProperties:
    @staticmethod
    def is_zero(matrix):
        return all(all(v == 0 for v in row) for row in matrix.rows)

    @staticmethod
    def is_symmetric(matrix):
        return MatrixProperties.is_square(matrix) and matrix.rows == matrix.transpose().rows

    @staticmethod
    def is_hermitian(matrix):
        return matrix.field == "complex" and matrix.rows == matrix.conjugate_transpose().rows

    @staticmethod
    def is_square(matrix):
        return len(matrix.rows) == len(matrix.rows[0])

    @staticmethod
    def determinant(matrix):
        if not MatrixProperties.is_square(matrix):
            raise ValueError("Determinant is only defined for square matrices.")
        if len(matrix.rows) == 1:
            return matrix.rows[0][0]
        return sum(((-1) ** col) * matrix.rows[0][col] *
                   MatrixProperties.determinant(Matrix(matrix.field,
                                                        [row[:col] + row[col + 1:] for row in matrix.rows[1:]]))
                   for col in range(len(matrix.rows)))


if __name__ == "__main__":
    m1 = Matrix("real", [[1, 0], [0, 1]])
    print("Matrix M1:")
    print(m1)
    print("\nIs M1 zero?", MatrixProperties.is_zero(m1))
    print("Is M1 symmetric?", MatrixProperties.is_symmetric(m1))
    print("Is M1 square?", MatrixProperties.is_square(m1))
    print("Determinant of M1:", MatrixProperties.determinant(m1))

    m2 = Matrix("real", [[2, 1], [1, 2]])
    print("\nMatrix M2:")
    print(m2)
    print("Is M2 symmetric?", MatrixProperties.is_symmetric(m2))
    print("Determinant of M2:", MatrixProperties.determinant(m2))
