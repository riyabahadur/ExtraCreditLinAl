class ComplexNumber:
    def __init__(self, re, im):
        self.real = re
        self.imag = im

    def __add__(self, other):
        return ComplexNumber(self.real + other.real, self.imag + other.imag)

    def __sub__(self, other):
        return ComplexNumber(self.real - other.real, self.imag - other.imag)

    def __mul__(self, other):
        re = self.real * other.real - self.imag * other.imag
        im = self.real * other.imag + self.imag * other.real
        return ComplexNumber(re, im)

    def conjugate(self):
        return ComplexNumber(self.real, -self.imag)

    def __eq__(self, other):
        return self.real == other.real and self.imag == other.imag

    def __str__(self):
        return f"{self.real} + {self.imag}i" if self.imag >= 0 else f"{self.real} - {-self.imag}i"


class Matrix:
    def __init__(self, fld, rws=None):
        self.field = fld
        self.rows = rws
        if rws and not all(len(r) == len(rws[0]) for r in rws):
            raise ValueError("All rows must have the same length.")

    def __mul__(self, other):
        if len(self.rows[0]) != len(other.rows):
            raise ValueError("Matrix dimensions do not align for multiplication.")
        res = [[sum(self.rows[i][k] * other.rows[k][j] for k in range(len(other.rows)))
                for j in range(len(other.rows[0]))] for i in range(len(self.rows))]
        return Matrix(self.field, rows=res)

    def transpose(self):
        transposed = [[self.rows[j][i] for j in range(len(self.rows))] for i in range(len(self.rows[0]))]
        return Matrix(self.field, rows=transposed)

    def conjugate(self):
        # Return the conjugate of the matrix (for complex fields only)
        if self.field != "complex":
            raise ValueError("Conjugate is only applicable for complex matrices.")
        conj = [[v.conjugate() if isinstance(v, ComplexNumber) else v for v in r] for r in self.rows]
        return Matrix(self.field, rows=conj)

    def conjugate_transpose(self):
        # Return the conjugate transpose of the matrix
        return self.conjugate().transpose()

    def __str__(self):
        # String representation of matrix
        return "\n".join([str(r) for r in self.rows])


class MatrixProperties:
    @staticmethod
    def is_zero(mat):
        return all(all(val == 0 for val in r) for r in mat.rows)

    @staticmethod
    def is_symmetric(mat):
        if not MatrixProperties.is_square(mat):
            return False
        return mat.rows == mat.transpose().rows

    @staticmethod
    def is_hermitian(mat):
        if mat.field != "complex":
            raise ValueError("Hermitian property is only applicable to complex matrices.")
        return mat.rows == mat.conjugate_transpose().rows

    @staticmethod
    def is_square(mat):
        return len(mat.rows) == len(mat.rows[0])

    @staticmethod
    def is_orthogonal(mat):
        if not MatrixProperties.is_square(mat):
            return False
        ident = Matrix("real", rows=[[1 if i == j else 0 for j in range(len(mat.rows))] for i in range(len(mat.rows))])
        return (mat * mat.transpose()).rows == ident.rows

    @staticmethod
    def is_unitary(mat):
        # Checking if the matrix is unitary (complex matrices only)
        if mat.field != "complex":
            raise ValueError("Unitary property is only applicable to complex matrices.")
        if not MatrixProperties.is_square(mat):
            return False
        ident = Matrix("complex", rows=[[1 if i == j else 0 for j in range(len(mat.rows))] for i in range(len(mat.rows))])
        return (mat * mat.conjugate_transpose()).rows == ident.rows

    @staticmethod
    def determinant(mat):
        # Calculate the determinant of the matrix (square matrices only)
        if not MatrixProperties.is_square(mat):
            raise ValueError("Determinant is only defined for square matrices.")
        if len(mat.rows) == 1:
            return mat.rows[0][0]
        det = 0
        for c in range(len(mat.rows)):
            sub_mat = Matrix("real", rows=[r[:c] + r[c+1:] for r in mat.rows[1:]])
            det += ((-1) ** c) * mat.rows[0][c] * MatrixProperties.determinant(sub_mat)
        return det


if __name__ == "__main__":
    m1 = Matrix("real", rows=[[1, 0], [0, 1]])
    print("Matrix M1:")
    print(m1)
    print("\nIs M1 zero?", MatrixProperties.is_zero(m1))
    print("Is M1 symmetric?", MatrixProperties.is_symmetric(m1))
    print("Is M1 square?", MatrixProperties.is_square(m1))
    print("Determinant of M1:", MatrixProperties.determinant(m1))

    m2 = Matrix("real", rows=[[2, 1], [1, 2]])
    print("\nMatrix M2:")
    print(m2)
    print("Is M2 symmetric?", MatrixProperties.is_symmetric(m2))
    print("Determinant of M2:", MatrixProperties.determinant(m2))
