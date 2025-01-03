import unittest
import ntru
import numpy as np
from polynomial_helper import Polynomial
class TestNtru(unittest.TestCase):
    def test_center(self):
        coefficients = [3,7,12]
        q = 16

        self.assertEqual([3,7,-4],ntru.center_coefficients(coefficients,q))

        coefficients = [-8,4,8,-5]
        q = 16

        self.assertEqual([-8,4,-8,-5],ntru.center_coefficients(coefficients,q))

    def test_encrypt_decrypt(self):
        N = 251
        p = 3
        q = 128
        d = 3

        for i in range(5):
            ntru_system = ntru.NTRU(N, p, q, d)
            public_key, private_key = ntru_system.generate_keys()
            message = Polynomial(np.random.choice([-1,0,1],N-1))
            ciphertext = ntru_system.encrypt(public_key, message)
            decrypted_message = ntru_system.decrypt(ciphertext, private_key)
            self.assertEqual(message, decrypted_message)

    def test_encode_message(self):
        # TODO: Implement
        pass