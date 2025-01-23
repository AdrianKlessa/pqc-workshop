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
        N = 509
        p = 3
        q = 2048
        d = 11

        for i in range(5):
            ntru_system = ntru.NTRU(N, p, q, d)
            public_key, private_key = ntru_system.generate_keys()
            message = Polynomial(np.random.choice([-1,0,1],N-1))
            ciphertext = ntru_system.encrypt(public_key, message)
            decrypted_message = ntru_system.decrypt(ciphertext, private_key)
            self.assertEqual(message, decrypted_message)

    def test_encode_decode_message(self):
        m1 = 0
        m2 = 127
        m3 = 8
        m4 = 7

        m1_encoded = ntru.encode_message(m1)
        m2_encoded = ntru.encode_message(m2)
        m3_encoded = ntru.encode_message(m3)
        m4_encoded = ntru.encode_message(m4)

        m1_decoded = ntru.decode_message(m1_encoded)
        m2_decoded = ntru.decode_message(m2_encoded)
        m3_decoded = ntru.decode_message(m3_encoded)
        m4_decoded = ntru.decode_message(m4_encoded)

        self.assertEqual(m1, m1_decoded)
        self.assertEqual(m2, m2_decoded)
        self.assertEqual(m3, m3_decoded)
        self.assertEqual(m4, m4_decoded)