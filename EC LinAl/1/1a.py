class ComplexNumber:
    def __init__(self, real, imag=0):
        self.real = real
        self.imag = imag

    def __add__(self, other):
        new_real = self.real + other.real
        new_imag = self.imag + other.imag
        return ComplexNumber(new_real, new_imag)

    def __mul__(self, other):
        real_part = (self.real * other.real) - (self.imag * other.imag)
        imag_part = (self.real * other.imag) + (self.imag * other.real)
        return ComplexNumber(real_part, imag_part)

    def __truediv__(self, other):
        if other.real == 0 and other.imag == 0:#divisibily by 0
            raise ZeroDivisionError("Cannot divide by zero.")
        
        denominator = (other.real**2) + (other.imag**2)
        
        # calulating real and imaginary parts
        real_part = ((self.real * other.real) + (self.imag * other.imag)) / denominator
        imag_part = ((self.imag * other.real) - (self.real * other.imag)) / denominator
        
        return ComplexNumber(real_part, imag_part)

    def abs(self):
        # to gwt the  magnitude of the complex number
        return (self.real**2 + self.imag**2)**0.5

    def conjugate(self):
        return ComplexNumber(self.real, -self.imag)

    def __str__(self):
        if self.imag >= 0:
            return f"{self.real} + {self.imag}i"
        else:
            return f"{self.real} - {-self.imag}i"
