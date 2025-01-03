from typing import Sequence, Tuple, Union
import numpy as np
import numpy.typing as npt
from polynomial_helper import multiplicative_inverse, Polynomial


# https://stackoverflow.com/a/2426520
def ntru_division(a: Polynomial, b: Polynomial, modulus: int) -> Tuple[
    Polynomial, Polynomial]:
    r = a
    q = Polynomial([0])
    u = multiplicative_inverse(int(b.coefficients[-1]), modulus)
    N = b.degree
    while r.degree >= N:
        d = r.degree
        x_coefficients = np.zeros(d+1)
        x_coefficients[d - N] = 1
        v = Polynomial([u]).multiply_mod(Polynomial([r.coefficients[d]])).multiply_mod(Polynomial(x_coefficients))
        v = v.reduced_modulo_scalar(modulus)
        r = r.substract_mod(v.multiply_mod(b))
        r = r.reduced_modulo_scalar(modulus)
        q = q.add_mod(v)
        q = q.reduced_modulo_scalar(modulus)
    return q, r


def extended_euclid(a: Polynomial, b: Polynomial, modulus: int) -> Tuple[Polynomial, Polynomial, Polynomial]:
    if b.is_zero:
        return Polynomial([1]), Polynomial([0]), a
    u = Polynomial([1])
    d = a
    v1 = Polynomial([0])
    v3 = b
    while not v3.is_zero:

        q, t3 = ntru_division(d, v3, modulus)
        t1 = u.substract_mod(q.multiply_mod(v1), modulus)
        u = v1
        d = v3
        v1 = t1
        v3 = t3
    to_divide = d.substract_mod(a.multiply_mod(u, modulus), modulus)
    v, r = ntru_division(to_divide, b, modulus)
    assert r.is_zero
    return (u, v, d)

def inverse_z_prime(a: Polynomial, modulus: int, modulus_polynomial: Polynomial) -> Union[Polynomial, False]:
    u, v, d = extended_euclid(a, modulus_polynomial, modulus)
    if d.degree==0:
        d_number = d.coefficients[0]
        d_inverse = multiplicative_inverse(d_number, modulus)
        b = Polynomial([d_inverse]).multiply_mod(u)
        return b
    return False

def inverse_z_prime_power(a: Polynomial, modulus_prime: int, modulus_power: int, modulus_polynomial: Polynomial):
    b = inverse_z_prime(a, modulus_prime, modulus_polynomial)
    if not b:
        return False
    n = modulus_prime
    while n<=modulus_power:
        modulus_number = modulus_prime**modulus_power
        b_squared = b.multiply_mod(b,modulus_number, modulus_polynomial)
        left = Polynomial([modulus_prime]).multiply_mod(b,modulus_number, modulus_polynomial)
        right = a.multiply_mod(b_squared,modulus_number, modulus_polynomial)
        b = left.substract_mod(right, modulus_number)
        n = n+1 #modulus_power*n
    _, r = ntru_division(b, modulus_polynomial, modulus_power)
    return r

def get_inverse(a: Polynomial, modulus_number: int, modulus_polynomial: Polynomial) -> Union[Polynomial, False]:
    if check_if_power_of_two(modulus_number):
        base = 2
        exponent = get_power_of_two(modulus_number)
        return inverse_z_prime_power(a, base, exponent, modulus_polynomial)
    else:
        return inverse_z_prime(a, modulus_number, modulus_polynomial)

def get_power_of_two(number: int) -> int:
    # int--> int
    # - 1 because 2ⁿ requires n+1 bits. Works for very large integers
    return number.bit_length()-1

def check_if_power_of_two(number: int) -> bool:
    return (number != 0) and (number & (number - 1) == 0)
