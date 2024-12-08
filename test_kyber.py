import unittest
import numpy as np
from polynomial_helper import Polynomial
import kyber
from kyber import KyberSystem


def all_elements_instance_of(array, cls):
    # Use np.vectorize to apply isinstance to each element
    check_isinstance = np.vectorize(lambda x: isinstance(x, cls))
    # Check all elements
    return np.all(check_isinstance(array))


class TestKyber(unittest.TestCase):
    def test_key_generation(self):
        ks = KyberSystem(256, 2, 3329, 2, 2)
        (A, t), s = ks.generate_key_pair()

        # Check the shape of arrays
        self.assertEqual(A.shape, (2, 2))
        self.assertEqual(t.shape, (2, 1))

        # Check that the elements are polynomials
        self.assertTrue(all_elements_instance_of(A, Polynomial))
        self.assertTrue(all_elements_instance_of(t, Polynomial))

    def test_encryption_decryption(self):
        ks = KyberSystem(256, 2, 3329, 2, 2)
        (A, t), s = ks.generate_key_pair()

        messages_to_encrypt = [1, 2, 3, 4, 5, 6, 7, 20, 100, 240, 256, 512, 520, 4920, 5192510, 5192511,
                               59113]  # Wide range of numbers
        for message in messages_to_encrypt:
            for i in range(10):  # Check for non-deterministic failures which should be incredibly rare
                u, v = ks.encrypt(A, t, message)
                decrypted = ks.decrypt(s, u, v)
                self.assertEqual(message, decrypted)

    def test_decimal_to_binary(self):
        num1 = 7
        num2 = 10
        num3 = 2
        self.assertEqual("111", kyber.decimal_to_binary(num1))
        self.assertEqual("1010", kyber.decimal_to_binary(num2))
        self.assertEqual("10", kyber.decimal_to_binary(num3))

    def test_get_nearest_to_half(self):
        q = 3329
        self.assertEqual(1665, kyber.get_nearest_integer_to_half(q))  # Python's round would return 1664

    def test_denoise(self):
        ks = KyberSystem(7, 2, 17, 2, 2)
        noisy_coefficients = [7, 14, 7, 5]
        actual_denoised = ks.denoise(noisy_coefficients)
        expected_denoised = [1, 0, 1, 1]
        self.assertEqual(expected_denoised, actual_denoised)

    def test_encode_message_to_polynomial(self):
        ks = KyberSystem(7, 2, 17, 2, 2)
        message_1 = 5  # 0b101
        message_2 = 8  # 0b1000

        encoded_1 = ks.encode_message_to_polynomial(message_1)
        expected_polynomial_1 = Polynomial([9, 0, 9])

        encoded_2 = ks.encode_message_to_polynomial(message_2)
        expected_polynomial_2 = Polynomial([0, 0, 0, 9])

        self.assertEqual(expected_polynomial_1, encoded_1)
        self.assertEqual(expected_polynomial_2, encoded_2)
