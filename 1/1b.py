class Vector:
    # Initializes a vector with a specified field (real or complex) and values
    def __init__(self, field, values):
        self.field = field
        self.values = values
        if not all(isinstance(v, (int, float, ComplexNumber)) for v in values):
            raise TypeError("Vector values must be real or complex numbers.")

    # Overloads the addition operator for vectors
    def __add__(self, other):
        if len(self.values) != len(other.values):
            raise ValueError("Vectors must have the same length for addition.")
        return Vector(self.field, [a + b for a, b in zip(self.values, other.values)])

    # Overloads scalar multiplication for vectors
    def __mul__(self, scalar):
        return Vector(self.field, [scalar * v for v in self.values])

    # Formats the vector as a string
    def __str__(self):
        return f"Vector({self.values})"

