from typing import Tuple

import numpy as np

import polynomial_helper
from polynomial_helper import Polynomial
from ntru_polynomials_wrapper import get_inverse


def get_ternary_polynomial(no_positive, no_negative, degree):
    """
    :param no_positive: Number of 1s in the vector
    :param no_negative: Number of -1s in the vector
    :param degree: Maximum degree of the vector
    :return:
    """
    vec = np.zeros(degree + 1)
    indices_positive = np.random.choice(degree + 1, no_positive, replace=False)
    vec[indices_positive] = 1
    possible_negative_indices = np.arange(degree + 1)
    possible_negative_indices = np.delete(possible_negative_indices, indices_positive)
    indices_negative = np.random.choice(possible_negative_indices, no_negative, replace=False)
    vec[indices_negative] = -1
    return Polynomial(vec)


def center_coefficients(coefficients, q):
    """
    Center the coefficients of a polynomial with respect to q, such that for each coefficient a
    q/2 <= a < q
    :param coefficients: Coefficients to center
    :param q: Center value
    :return: Coefficients centered to range [-q/2,q/2)
    """
    centered = []
    for c in coefficients:
        if c >= q / 2:
            centered.append((-q / 2) + (c % (q / 2)))
        elif c < -q / 2:
            centered.append((c % (q / 2)))
        else:
            centered.append(c)
    return centered

def center_polynomial(polynomial: Polynomial, q: int):
    coeffs = polynomial.coefficients
    return Polynomial(center_coefficients(coeffs, q))

class NTRU:
    def __init__(self, N: int = 251, p: int = 3, q: int = 128, d: int = 3):
        self.N = N
        self.p = p
        self.q = q
        self.d = d
        self.N_polynomial = polynomial_helper.Polynomial([-1] + [0] * (N - 1) + [1])

    def generate_keys(self):
        while True:
            # Repeat until invertible
            # Some discussion about the choice of parameters: https://crypto.stackexchange.com/a/2634
            # so for ternary(d,d) the chance of it being invertible might be low
            f = get_ternary_polynomial(self.d, self.d-1, self.N - 1)
            if f.evaluated_at_1 == 0:
                continue
            fp = get_inverse(f, self.p, self.N_polynomial)
            fq = get_inverse(f, self.q, self.N_polynomial)
            if fp and fq:
                break
        g = get_ternary_polynomial(self.d, self.d, self.N - 1)
        h = Polynomial([self.p]).multiply_mod(fq, self.q, self.N_polynomial)
        h = h.multiply_mod(g, self.q, self.N_polynomial)
        h = center_polynomial(h, self.q)
        return (h, (f, fp))

    def encrypt(self, public_key: Polynomial, message: Polynomial) -> Polynomial:
        for coeff in message.coefficients:
            if not -self.p/2 <= coeff <= self.p/2:
                raise ValueError("Message is not centered properly")
        r = get_ternary_polynomial(self.d, self.d, self.N - 1)
        c = r.multiply_mod(public_key, self.q, self.N_polynomial)
        c = c.add_mod(message, self.q)
        return c

    def decrypt(self, ciphertext: Polynomial, private_key: Tuple[Polynomial, Polynomial]) -> Polynomial:
        a = ciphertext.multiply_mod(private_key[0], self.q, self.N_polynomial)
        a = center_polynomial(a, self.q)
        m = private_key[1].multiply_mod(a, self.p, self.N_polynomial)
        m = center_polynomial(m, self.p)
        return m

def encode_message(message: int)->Polynomial:
    binary = str(bin(message))[2:]
    binary = [int(i) for i in binary]
    # Reversing the binary strings to prevent truncation of leading coefficients equal to 0
    binary = binary[::-1]
    return Polynomial(binary)

def decode_message(message_polynomial: Polynomial)-> int:
    coefficients = message_polynomial.coefficients
    binary = coefficients[::-1].tolist()
    binary = [str(i) for i in binary]
    binary = ''.join(binary)
    return int(binary, 2)