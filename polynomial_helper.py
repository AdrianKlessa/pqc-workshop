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
            # Trim zero coefficients to avoid equality, other issues
            self.coefficients = np.trim_zeros(self.coefficients, 'b')

    def multiply_mod(self, other_polynomial: Polynomial, mod_number: int | None = None,
                     mod_polynomial: Polynomial | None = None) -> Polynomial:
        """
        Multiply this by other polynomial. Optionally reduce modulo.
        :param other_polynomial: Polynomial to multiply by
        :param mod_number: Modulus for coefficients
        :param mod_polynomial: Modulus polynomial for reduction
        :return: this*other_polynomial
        """
        result_coefficients = np.convolve(self.coefficients, other_polynomial.coefficients)
        if mod_polynomial is not None:
            _, result = Polynomial(result_coefficients).divide_by(mod_polynomial, mod_number)
            return result
        if mod_number is not None:
            result_coefficients %= mod_number
            return Polynomial(result_coefficients)

        return Polynomial(result_coefficients)

    def add_mod(self, other_polynomial: Polynomial, mod_number: int | None = None) -> Polynomial:
        """
        Add another polynomial to this polynomial. Optionally reduce modulo.
        :param other_polynomial: Polynomial to add to
        :param mod_number: Modulus for coefficients
        :return: this+other_polynomial
        """
        # this + other
        len_a = len(self.coefficients)
        len_b = len(other_polynomial.coefficients)
        if len_a > len_b:
            result = self.coefficients.copy()
            result[:len_b] += other_polynomial.coefficients
        else:
            result = other_polynomial.coefficients.copy()
            result[:len_a] += self.coefficients

        if mod_number:
            result %= mod_number
        return Polynomial(result)

    def substract_mod(self, other_polynomial: Polynomial, mod_number: int | None = None) -> Polynomial:
        """
        Subtract another polynomial from this polynomial. Optionally reduce modulo.
        :param other_polynomial: Polynomial to subtract from this polynomial
        :param mod_number: Modulus for coefficients
        :return: this-other_polynomial
        """
        temp_coeffs = np.copy(other_polynomial.coefficients)
        temp_coeffs *= -1
        return self.add_mod(Polynomial(temp_coeffs), mod_number)

    def reduced_modulo_scalar(self, scalar: int) -> Polynomial:
        """
        Reduce this polynomial's coefficients modulo scalar
        :param scalar: Modulus for coefficients
        :return: Polynomial with reduced coefficients
        """
        coeffs = np.copy(self.coefficients)
        coeffs = coeffs % scalar
        return Polynomial(coeffs)

    def divide_by(self, other_polynomial: Polynomial, mod_number: int | None = None) -> tuple[Polynomial, Polynomial]:
        """
        Divide this polynomial by another polynomial. Optionally reduce modulo.
        :param other_polynomial: Polynomial to divide by
        :param mod_number: Modulus for coefficients
        :return: Tuple of (quotient, remainder)
        """
        this = self.reduced_modulo_scalar(mod_number) if mod_number else self
        other = other_polynomial.reduced_modulo_scalar(mod_number) if mod_number else other_polynomial
        if other.is_zero:
            raise ValueError("Cannot divide by zero")
        if this == other:
            return Polynomial([1]), Polynomial([0])
        # https://en.wikipedia.org/wiki/Polynomial_greatest_common_divisor#Euclidean_division
        q = Polynomial([0])
        r = this
        d = other.degree
        c = other.coefficients[-1]
        if this.degree > other.degree:
            deg_var = r.degree
            while deg_var >= d and not r.is_zero:
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
        else:
            q = Polynomial([0])
            r = this
        if mod_number:
            q = q.reduced_modulo_scalar(mod_number)
            r = r.reduced_modulo_scalar(mod_number)

        assert other_polynomial.multiply_mod(q, mod_number).add_mod(r, mod_number) == this
        return q, r

    def get_inverse(self, mod_polynomial: Polynomial, mod_number: int | None) -> Polynomial:

        """
        Get the multiplicative inverse of this polynomial in a ring modulo mod_polynomial
        :param mod_polynomial: Modulus polynomial
        :param mod_number: Modulus for coefficients
        :return: Multiplicative inverse of this polynomial
        """
        this = self.reduced_modulo_scalar(mod_number)
        other = mod_polynomial.reduced_modulo_scalar(mod_number)
        # https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm#Simple_algebraic_field_extensions
        t = Polynomial([0])
        t_prime = Polynomial([1])
        r = other
        r_prime = this
        while not r_prime.is_zero:
            q, _ = r.divide_by(r_prime, mod_number)
            left = r_prime
            right = r.substract_mod(q.multiply_mod(r_prime), mod_number)
            r = left
            r_prime = right
            left = t_prime
            right = t.substract_mod(q.multiply_mod(t_prime), mod_number)
            t = left
            t_prime = right
        _, r = r.divide_by(mod_polynomial, mod_number)
        if r.degree > 0:
            raise ValueError("The given polynomial is not invertible")
        temp_x = multiplicative_inverse(int(r.coefficients[0]), mod_number)
        return Polynomial([temp_x]).multiply_mod(t, mod_number)

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

    @property
    def evaluated_at_1(self):
        return np.sum(self.coefficients)

    def __eq__(self, other):
        return np.array_equal(self.coefficients, other.coefficients)

    def __str__(self):
        return repr(self)

    # Algebraic operation overloads for plug & play Numpy matrix operations

    def __add__(self, other):
        return self.add_mod(other)

    def __sub__(self, other):
        return self.substract_mod(other)

    def __mul__(self, other):
        return self.multiply_mod(other)


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
