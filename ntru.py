from typing import Tuple

import numpy as np

import polynomial_helper
from polynomial_helper import Polynomial


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


class NTRU:
    def __init__(self, N:int, p:int, q:int, d: int):
        self.N = N
        self.p = p
        self.q = q
        self.d = d
        self.N_polynomial = polynomial_helper.Polynomial([-1]+[0]*(N-1)+[1])
        print(self.N_polynomial)

    def generate_keys(self):
        while True:
            try:
                # Repeat until invertible
                print("Aaa")
                f = get_ternary_polynomial(self.d, self.d, self.N - 1)
                fp = f.get_inverse(self.N_polynomial, self.p)
                print("Got fp")
                fq = f.get_inverse(self.N_polynomial, self.q)
                break
            except ValueError:
                print("retrying...")
                continue
        g = get_ternary_polynomial(self.d,self.d,self.N-1)
        h = Polynomial([self.p]).multiply_mod(fq,self.q,self.N_polynomial)
        h = h.multiply_mod(g,self.q,self.N_polynomial)
        return (h,(f, fp))

    def encrypt(self, public_key: Polynomial, message: Polynomial)->Polynomial:
        r = get_ternary_polynomial(self.d, self.d, self.N - 1)
        c = r.multiply_mod(public_key, self.q, self.N_polynomial)
        c = c.multiply_mod(Polynomial[self.p], self.q, self.N_polynomial)
        c = c.add_mod(message, self.q)
        return c

    def decrypt(self, ciphertext: Polynomial, private_key: Tuple[Polynomial,Polynomial])->Polynomial:
        a = ciphertext.multiply_mod(private_key[0], self.q, self.N_polynomial)
        m = private_key[1].multiply_mod(a, self.p, self.N_polynomial)
        return m


ntru_system = NTRU(11,3,32, 3)

public_key, private_key = ntru_system.generate_keys()

message = Polynomial([-1,0,0,1,-1,0,0,0,-1,1,1])
print(f"message: {message}")

ciphertext = ntru_system.encrypt(public_key, message)

print(f"ciphertext: {ciphertext}")

decrypted_message = ntru_system.decrypt(ciphertext, private_key)

print(f"decrypted_message: {decrypted_message}")