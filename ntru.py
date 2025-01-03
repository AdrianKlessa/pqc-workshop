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
    def __init__(self, N: int, p: int, q: int, d: int):
        self.N = N
        self.p = p
        self.q = q
        self.d = d
        self.N_polynomial = polynomial_helper.Polynomial([-1] + [0] * (N - 1) + [1])
        print(self.N_polynomial)

    def generate_keys(self):
        while True:
            # Repeat until invertible
            # Some discussion about the choice of parameters: https://crypto.stackexchange.com/a/2634
            # so for ternary(d,d) the chance of it being might be low
            f = get_ternary_polynomial(self.d, self.d-1, self.N - 1)
            if f.evaluated_at_1 == 0:
                continue
            f = Polynomial([-1, 1, 1, 0, -1, 0, 1, 0, 0, 1, -1]) # TODO: Remove hard-coded value
            print(f"f: {f}")
            print(f"N polynomial: {self.N_polynomial}")
            print(f"q: {self.q}")
            print(f"p: {self.p}")
            fp = get_inverse(f, self.p, self.N_polynomial)
            fq = get_inverse(f, self.q, self.N_polynomial)
            print(f"fp: {fp}, fq: {fq}")
            if fp and fq:
                break
        # TODO: Remove hard-coded value
        #g = get_ternary_polynomial(self.d, self.d, self.N - 1)
        g = Polynomial([-1,0,1,1,0,1,0,0,-1,0,-1])
        print(f"g: {g}")
        h = Polynomial([self.p]).multiply_mod(fq, self.q, self.N_polynomial)
        h = h.multiply_mod(g, self.q, self.N_polynomial)
        h = center_polynomial(h, self.q)
        print(f"h: {h}")
        return (h, (f, fp))

    def encrypt(self, public_key: Polynomial, message: Polynomial) -> Polynomial:
        # TODO: Verify that the message is in the correct range (centered)
        # TODO: Remove hard-coded value
        #r = get_ternary_polynomial(self.d, self.d, self.N - 1)
        r = Polynomial([-1,0,1,1,1,-1,0,-1])
        print(f"r: {r}")
        print(f"public_key: {public_key}")
        c = r.multiply_mod(public_key, self.q, self.N_polynomial)
        #c = c.multiply_mod(Polynomial([self.p]), self.q, self.N_polynomial)
        c = c.add_mod(message, self.q)
        return c

    def decrypt(self, ciphertext: Polynomial, private_key: Tuple[Polynomial, Polynomial]) -> Polynomial:
        a = ciphertext.multiply_mod(private_key[0], self.q, self.N_polynomial)
        a = center_polynomial(a, self.q)
        m = private_key[1].multiply_mod(a, self.p, self.N_polynomial)
        m = center_polynomial(m, self.p)
        return m

# TODO: Force a specific f polynomial and use the one from wikipedia for testing as it's small
#ntru_system = NTRU(509, 3, 2048, 3)
ntru_system = NTRU(11, 3, 32, 3)
public_key, private_key = ntru_system.generate_keys()

message = Polynomial([-1, 0, 0, 1, -1, 0, 0, 0, -1, 1, 1])
print(f"message: {message}")

ciphertext = ntru_system.encrypt(public_key, message)

print(f"ciphertext: {ciphertext}")

decrypted_message = ntru_system.decrypt(ciphertext, private_key)

print(f"decrypted_message: {decrypted_message}")

#assert decrypted_message.reduced_modulo_scalar(3) == message.reduced_modulo_scalar(3)

assert decrypted_message == message