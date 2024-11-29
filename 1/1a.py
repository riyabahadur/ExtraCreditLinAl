class ComplexNumber:
    # Initializes a complex number with real and imaginary parts
    def __init__(self, real, imag=0):
        self.real = real
        self.imag = imag

    # Overloads the addition operator for complex numbers
    def __add__(self, other):
        return ComplexNumber(self.real + other.real, self.imag + other.imag)

    # Overloads the multiplication operator for complex numbers
    def __mul__(self, other):
        real_part = self.real * other.real - self.imag * other.imag
        imag_part = self.real * other.imag + self.imag * other.real
        return ComplexNumber(real_part, imag_part)

    # Overloads the division operator for complex numbers
    def __truediv__(self, other):
        if other.real == 0 and other.imag == 0:
            raise ZeroDivisionError("Cannot divide by zero.")
        denominator = other.real**2 + other.imag**2
        real_part = (self.real * other.real + self.imag * other.imag) / denominator
        imag_part = (self.imag * other.real - self.real * other.imag) / denominator
        return ComplexNumber(real_part, imag_part)

    # Returns the absolute value of the complex number
    def abs(self):
        return (self.real**2 + self.imag**2)**0.5

    # Returns the complex conjugate
    def conjugate(self):
        return ComplexNumber(self.real, -self.imag)


    def __str__(self):
        return f"{self.real} + {self.imag}i"
