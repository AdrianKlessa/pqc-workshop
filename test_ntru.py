import unittest
import ntru

class TestNtru(unittest.TestCase):
    def test_center(self):
        coefficients = [3,7,12]
        q = 16

        self.assertEqual([3,7,-4],ntru.center_coefficients(coefficients,q))

        coefficients = [-8,4,8,-5]
        q = 16

        self.assertEqual([-8,4,-8,-5],ntru.center_coefficients(coefficients,q))