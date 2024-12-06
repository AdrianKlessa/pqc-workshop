from __future__ import annotations

from typing import Sequence

import numpy as np


class Polynomial:

    def __init__(self, coefficients: Sequence[int]):
        """

        :param coefficients: Coefficients, starting from the lowest degree (a0, a1, a2, ...)
        """
        if not np.any(coefficients):
            self.coefficients = np.array([0], dtype=int)
        else:
            self.coefficients = np.array(coefficients, dtype=int)
            self.coefficients = np.trim_zeros(self.coefficients, 'b')

    def multiply_mod(self, other_polynomial: Polynomial, mod_number: int | None = None,
                     mod_polynomial: Polynomial | None = None) -> Polynomial:
        # +2 since the degree is w/o the constant term
        result_coefficients = [0 for _ in
                               range(self.coefficients.shape[0] + other_polynomial.coefficients.shape[0] + 2)]

        for i in range(other_polynomial.coefficients.shape[0]):
            for j in range(self.coefficients.shape[0]):
                result_coefficients[i + j] += (self.coefficients[j] * other_polynomial.coefficients[i])
                if mod_number:
                    result_coefficients[i + j] %= mod_number

        if np.any(result_coefficients):
            last_nonzero_index = np.max(np.nonzero(result_coefficients))
            result_coefficients = result_coefficients[:last_nonzero_index + 1]
        else:
            result_coefficients = np.array([0], dtype=int)

        if mod_polynomial:
            temp_polynomial = Polynomial(result_coefficients)
            q, r = temp_polynomial.divide_by(mod_polynomial, mod_number)
            return r

        return Polynomial(result_coefficients)

    def add_mod(self, other_polynomial: Polynomial, mod_number: int | None = None) -> Polynomial:
        # this + other
        if self.degree > other_polynomial.degree:
            coeffs = np.copy(self.coefficients)
            for i in range(other_polynomial.degree + 1):
                coeffs[i] = coeffs[i] + other_polynomial.coefficients[i]
        else:
            coeffs = np.copy(other_polynomial.coefficients)
            for i in range(self.degree + 1):
                coeffs[i] = coeffs[i] + self.coefficients[i]

        if np.any(coeffs):
            last_nonzero_index = np.max(np.nonzero(coeffs))
            coeffs = coeffs[:last_nonzero_index + 1]
        else:
            coeffs = np.array([0], dtype=int)

        if mod_number:
            for i in range(len(coeffs)):
                coeffs[i] %= mod_number
        return Polynomial(coeffs)

    def substract_mod(self, other_polynomial: Polynomial, mod_number: int | None = None) -> Polynomial:
        temp_coeffs = np.copy(other_polynomial.coefficients)
        temp_coeffs *= -1
        return self.add_mod(Polynomial(temp_coeffs), mod_number)

    def reduced_modulo_scalar(self, scalar: int) -> Polynomial:
        coeffs = np.copy(self.coefficients).tolist()
        for i in range(len(coeffs)):
            coeffs[i] %= scalar
        return Polynomial(coeffs)

    def divide_by(self, other_polynomial: Polynomial, mod_number: int | None = None) -> tuple[Polynomial, Polynomial]:
        # this / other
        # --> quotient, remainder
        """
        to_divide = self
        while to_divide.degree>=other_polynomial.degree:
            to_divide
        """

        if other_polynomial.is_zero:
            raise ValueError("Cannot divide by zero")

        this = self.reduced_modulo_scalar(mod_number) if mod_number else self
        other = other_polynomial.reduced_modulo_scalar(mod_number) if mod_number else other_polynomial
        # https://en.wikipedia.org/wiki/Polynomial_greatest_common_divisor#Euclidean_division
        q = Polynomial([0])
        r = this
        d = other.degree
        c = other.coefficients[-1]

        deg_var = r.degree
        while deg_var >= d:
            if mod_number:
                c_inverse = multiplicative_inverse(c, mod_number)
                s_left = (r.coefficients[-1] * c_inverse) % mod_number
            else:
                s_left = r.coefficients[-1] / c
            s_degree = r.degree - d
            s_coefficients = [0] * s_degree + [s_left]
            s = Polynomial(s_coefficients)
            q = q.add_mod(s, mod_number)
            r = r.substract_mod(s.multiply_mod(other, mod_number, None), mod_number)
            deg_var = r.degree
        if mod_number:
            q = q.reduced_modulo_scalar(mod_number)
            r = r.reduced_modulo_scalar(mod_number)
        return q, r

    def get_inverse(self, mod_polynomial: Polynomial, mod_number: int | None) -> Polynomial:

        this = self.reduced_modulo_scalar(mod_number) if mod_number else self
        _, this = this.divide_by(mod_polynomial, mod_number)
        other = mod_polynomial.reduced_modulo_scalar(mod_number) if mod_number else mod_polynomial


        # https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm#Simple_algebraic_field_extensions
        t = Polynomial([0])
        t_prime = Polynomial([1])
        r = mod_polynomial
        r_prime = self
        while not r_prime.is_zero:
            print(f"r: {r}")
            print(f"r_prime: {r_prime}")
            q, _ = r.divide_by(r_prime, mod_number)
            print(f"q: {q}")
            left = r_prime
            right = r.substract_mod(q.multiply_mod(r_prime, mod_number, None), mod_number)
            print(f"right1: {right}")
            print(f"left1: {left}")
            r = left
            r_prime = right
            left = t_prime
            right = t.substract_mod(q.multiply_mod(t_prime, mod_number, None), mod_number)
            print(f"right2: {right}")
            print(f"left2: {left}")
            t = left
            t_prime = right
            print(r_prime)
        if r.degree > 0:
            print(f"r degree: {r.degree}")
            print(f"r: {r}")
            raise ValueError("The given polynomial is not invertible")
        return r.divide_by(Polynomial([1]), mod_number)[0].multiply_mod(t, mod_number,None)

    @property
    def degree(self):
        if self.is_zero:
            return -1
        return self.coefficients.shape[0] - 1

    @property
    def is_zero(self):
        return len(self.coefficients) == 0 or np.all(self.coefficients == 0)

    def __repr__(self):

        if all(v == 0 for v in self.coefficients):
            return "0"

        repr_string = ""
        for i in reversed(range(self.degree + 1)):
            coefficient_i = self.coefficients[i]

            # I hate this
            if coefficient_i == 0:
                continue
            if i == self.degree and i > 1:
                repr_string += f"{coefficient_i}x^{i}"
            elif i == self.degree and i == 1:
                repr_string += f"{coefficient_i}x"
            elif i == self.degree and i == 0:
                repr_string += f"{coefficient_i}"
            else:
                if coefficient_i > 0:
                    repr_string += "+"
                if i == 0:
                    repr_string += f"{coefficient_i}"
                elif i == 1:
                    repr_string += f"{coefficient_i}x"
                else:
                    repr_string += f"{coefficient_i}x^{i}"

        return repr_string

    def __eq__(self, other):
        return np.array_equal(self.coefficients, other.coefficients)

    def __str__(self):
        return repr(self)


def extendedGCD(a, b):
    r, r1 = a, b
    s, s1 = 1, 0
    t, t1 = 0, 1
    while r1 != 0:
        q, r2 = r // r1, r % r1
        r, s, t, r1, s1, t1 = r1, s1, t1, r2, s - s1 * q, t - t1 * q
    d = r
    return d, s, t


def multiplicative_inverse(a: int, m: int) -> int:
    d, inv, _ = extendedGCD(a, m)
    if d == 1:
        if m == 1:
            return 1  #for compatibility
        return inv % m
    else:
        raise ValueError('Numbers ' + str(a) + ' and ' + str(m) + ' are not coprime.')

#print(multiplicative_inverse(3,7))