class ComplexNumber:
    def __init__(self, real=0, imag=0):
        self.real = real
        self.imag = imag

    def __add__(self, other):
        if isinstance(other, ComplexNumber):
            return ComplexNumber(self.real + other.real, self.imag + other.imag)
        return NotImplemented

    def __mul__(self, other):
        if isinstance(other, ComplexNumber):
            real_part = self.real * other.real - self.imag * other.imag
            imag_part = self.real * other.imag + self.imag * other.real
            return ComplexNumber(real_part, imag_part)
        return NotImplemented

    def __repr__(self):
        return f"{self.real} + {self.imag}i"


class Vector:
    def __init__(self, field, n):
        self.field = field
        self.values = []
        if field == 'complex':
            for i in range(n):
                real = float(input(f"Enter real part for value {i+1}: "))
                imag = float(input(f"Enter imaginary part for value {i+1}: "))
                self.values.append(ComplexNumber(real, imag))
        elif field == 'real':
            for i in range(n):
                value = float(input(f"Enter value for coordinate {i+1}: "))
                self.values.append(value)
        else:
            raise ValueError("Field must be 'real' or 'complex'.")

    def __add__(self, other):
        # Add two vectors element-wise if they have the same length

        if len(self.values) != len(other.values):
            raise ValueError("Vectors must have the same length for addition.")
        return Vector(self.field, [a + b for a, b in zip(self.values, other.values)])

    def __mul__(self, scalar):
        # Multiply all vector elements by a scalar
        return Vector(self.field, [scalar * v for v in self.values])

    def __str__(self):
        return f"Vector({self.values})"


vec1 = Vector('real', 3)

vec2 = Vector('complex', 3)

print(vec1)  # Output: Vector([real values input by user])
print(vec2)  # Output: Vector([complex values input by user])
