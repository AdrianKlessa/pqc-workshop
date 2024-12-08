import unittest
import numpy as np
from polynomial_helper import Polynomial
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
                u, v = ks.encrypt(A,t, message)
                decrypted = ks.decrypt(s, u, v)
                self.assertEqual(message,decrypted)