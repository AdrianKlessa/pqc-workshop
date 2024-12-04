from __future__ import annotations

from typing import Sequence

import numpy as np


class Polynomial:
    def __init__(self, coefficients: Sequence[int]):
        """

        :param coefficients: Coefficients, starting from the lowest degree (a0, a1, a2, ...)
        """
        self.coefficients = np.array(coefficients)

    def multiply_mod(self, other_polynomial: Polynomial, mod_number: int | None,
                     mod_polynomial: Polynomial | None) -> Polynomial:
        # +2 since the degree is w/o the constant term
        result_coefficients = [0 for _ in range(self.degree + other_polynomial.degree + 2)]

        for i in range(other_polynomial.coefficients.shape[0]):
            for j in range(self.coefficients.shape[0]):
                result_coefficients[i + j] += self.coefficients[j] * other_polynomial.coefficients[i]

        if all(v == 0 for v in result_coefficients):
            return Polynomial([0])
        last_nonzero_index = np.max(np.nonzero(result_coefficients))
        result_coefficients = result_coefficients[:last_nonzero_index + 1]
        # TODO: Reduce mod_number and mod_polynomial
        return Polynomial(result_coefficients)

    def add_mod(self, other_polynomial: Polynomial, mod_number: int | None) -> Polynomial:
        # this + other
        # TODO: Reduce mod number
        if self.degree > other_polynomial.degree:
            coeffs = np.copy(self.coefficients)
            for i in range(other_polynomial.degree + 1):
                coeffs[i] = coeffs[i] + other_polynomial.coefficients[i]
        else:
            coeffs = np.copy(other_polynomial.coefficients)
            for i in range(self.degree + 1):
                coeffs[i] = coeffs[i] + self.coefficients[i]

        if all(v == 0 for v in coeffs):
            return Polynomial([0])
        last_nonzero_index = np.max(np.nonzero(coeffs))
        coeffs = coeffs[:last_nonzero_index + 1]
        return Polynomial(coeffs)

    def substract_mod(self, other_polynomial: Polynomial, mod_number: int | None) -> Polynomial:
        temp_coeffs = np.copy(other_polynomial.coefficients)
        temp_coeffs *= -1
        return self.add_mod(Polynomial(temp_coeffs), mod_number)

    def divide_by(self, other_polynomial: Polynomial, mod_number: int | None) -> (Polynomial, Polynomial):
        # this / other
        """
        to_divide = self
        while to_divide.degree>=other_polynomial.degree:
            to_divide
        """
        # https://en.wikipedia.org/wiki/Polynomial_greatest_common_divisor#Euclidean_division
        q = Polynomial([0])
        r = self
        d = other_polynomial.degree
        c = other_polynomial.coefficients[-1]

        deg_var = r.degree
        while deg_var >= d:
            s_left = r.coefficients[-1] / c
            s_degree = r.degree - d
            s_coefficients = [0] * s_degree + [s_left]
            s = Polynomial(s_coefficients)
            q = q.add_mod(s, mod_number)
            r = r.substract_mod(s.multiply_mod(other_polynomial, mod_number, None), mod_number)
            deg_var = r.degree
            #print(f"s: {s}")
            #print(f"q: {q}")
            #print(f"r: {r}")
            #print(deg_var)
            #raise Exception
        return q, r

    @property
    def degree(self):
        return self.coefficients.shape[0] - 1

    def __repr__(self):

        repr_string = ""
        for i in reversed(range(self.degree + 1)):
            coeff = self.coefficients[i]

            # I hate this
            if coeff == 0:
                continue
            if i == self.degree and i > 1:
                repr_string += f"{coeff}x^{i}"
            elif i == self.degree and i == 1:
                repr_string += f"{coeff}x"
            elif i == self.degree and i == 0:
                repr_string += f"{coeff}"
            else:
                if coeff>0:
                    repr_string +="+"
                if i == 0:
                    repr_string += f"{coeff}"
                elif i == 1:
                    repr_string += f"{coeff}x"
                else:
                    repr_string += f"{coeff}x^{i}"

        return repr_string

    def __str__(self):
        return repr(self)
