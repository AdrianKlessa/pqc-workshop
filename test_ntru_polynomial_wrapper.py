import unittest
from polynomial_helper import Polynomial
from ntru_polynomials_wrapper import get_inverse

class TestNTRUPolynomial(unittest.TestCase):
    def test_inverse_ntru(self):
        f = Polynomial([1, 0, -1, 1])
        p = 3
        q = 16
        N_polynomial = Polynomial([-1, 0, 0, 0, 1])

        expected_p = Polynomial([2, 0, 1, 1]).reduced_modulo_scalar(p)
        expected_q = Polynomial([-3, 7, 3, -6]).reduced_modulo_scalar(q)
        # Sanity check
        assert expected_p.multiply_mod(f, p, N_polynomial) == Polynomial([1])
        assert expected_q.multiply_mod(f, q, N_polynomial) == Polynomial([1])

        #f_p = f.get_inverse(N_polynomial, p)
        #f_q = f.get_inverse(N_polynomial, q)
        f_p = get_inverse(f, p, N_polynomial)
        f_q = get_inverse(f, q, N_polynomial)

        self.assertEqual(f_p, expected_p)
        self.assertEqual(f_q, expected_q)

    def test_inverse_ntru_2(self):
        # Example from wikipedia for NTRUEncrypt
        f = Polynomial([-1, 1, 1, 0, -1, 0, 1, 0, 0, 1, -1])  # -1+x+x^2-x^4+x^6+x^9-x^10
        p = 3
        q = 32
        N = 11
        modulus = Polynomial([-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1])
        f_p = Polynomial([1, 2, 0, 2, 2, 1, 0, 2, 1, 2])
        f_q = Polynomial([5, 9, 6, 16, 4, 15, 16, 22, 20, 18, 30])
        # print(f"inverse p: {f.reduced_modulo_scalar(p).get_inverse(modulus.reduced_modulo_scalar(p), p)}")
        # print(f"inverse q: {f.reduced_modulo_scalar(q).get_inverse(modulus.reduced_modulo_scalar(q), q)}")
        # Sanity check
        assert f_p.multiply_mod(f, p, modulus) == Polynomial([1])
        assert f_q.multiply_mod(f, q, modulus) == Polynomial([1])
        #a1 = f.get_inverse(modulus, p)
        #a2 = f.get_inverse(modulus, q)
        a1 = get_inverse(f, p, modulus)
        a2 = get_inverse(f, q, modulus)
        self.assertEqual(a1, f_p)
        self.assertEqual(a2, f_q)