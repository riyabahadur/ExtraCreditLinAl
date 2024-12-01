class Matrix:
    def __init__(self, field, rows=None, columns=None):
        self.field = field
        if rows:
            self.rows = rows
            if not all(len(row) == len(rows[0]) for row in rows):#checking if the rows have the same length
                raise ValueError("All rows must have the same length.")
        elif columns:
            self.rows = [list(col) for col in zip(*columns)]#transposing colums to form rows
        else:
            raise ValueError("Either rows or columns must be provided to initialize the matrix.")

    def __add__(self, other):
        if len(self.rows) != len(other.rows) or len(self.rows[0]) != len(other.rows[0]):
            raise ValueError("Matrix dimensions do not match for addition.")
        result = [[self.rows[i][j] + other.rows[i][j] for j in range(len(self.rows[0]))] for i in range(len(self.rows))]
        return Matrix(self.field, rows=result)

    def __mul__(self, other):
        if len(self.rows[0]) != len(other.rows):
            raise ValueError("Matrix dimensions do not align for multiplication.")
        result = [[sum(self.rows[i][k] * other.rows[k][j] for k in range(len(other.rows)))
                   for j in range(len(other.rows[0]))] for i in range(len(self.rows))]
        return Matrix(self.field, rows=result)

    def transpose(self):
        #returin the transpose of the matrix
        return Matrix(self.field, rows=[[self.rows[j][i] for j in range(len(self.rows))] for i in range(len(self.rows[0]))])

    def conjugate(self):
        if self.field != "complex":
            raise ValueError("Conjugate is only applicable for complex matrices.")
        return Matrix(self.field, rows=[[v.conjugate() if isinstance(v, ComplexNumber) else v for v in row] for row in self.rows]) # type: ignore

    def conjugate_transpose(self):
        return self.conjugate().transpose()

    def get_row(self, index):
        if index < 0 or index >= len(self.rows):
            raise IndexError("Row index out of bounds.")
        return self.rows[index]

    def get_column(self, index):
        if index < 0 or index >= len(self.rows[0]):
            raise IndexError("Column index out of bounds.")
        return [self.rows[i][index] for i in range(len(self.rows))]

    def __str__(self):
        return "\n".join([str(row) for row in self.rows])
