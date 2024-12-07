import unittest
from polynomial_helper import Polynomial, multiplicative_inverse
import numpy as np

class TestPolyMatrix(unittest.TestCase):
    def test_matrix_matrix_multiplication(self):
        a1 = Polynomial([0,1])
        a2 = Polynomial([0,0,1])
        a3 = Polynomial([1])
        a4 = Polynomial([0,1])

        b1 = Polynomial([0,1])
        b2 = Polynomial([0, 1])
        b3 = Polynomial([0, 1])
        b4 = Polynomial([1])

        c1 = Polynomial([0,0,1,1])
        c2 = Polynomial([0,0,2])
        c3 = Polynomial([0,1,1])
        c4 = Polynomial([0,2])

        A = np.array([[a1,a2],[a3,a4]], dtype=object)
        B = np.array([[b1,b2],[b3,b4]], dtype=object)
        C = np.array([[c1,c2],[c3,c4]], dtype=object)

        result = np.dot(A,B)

        self.assertEqual(True, np.array_equal(C,result))
    def test_matrix_vector_multiplication(self):
        a1 = Polynomial([0,1])
        a2 = Polynomial([0,0,1])
        a3 = Polynomial([1])
        a4 = Polynomial([0,1])

        b1 = Polynomial([1])
        b2 = Polynomial([0,1])

        A = np.array([[a1,a2],[a3,a4]], dtype=object)
        x = np.array([b1,b2],ndmin=2).T

        c1 = Polynomial([0,1,0,1])
        c2 = Polynomial([1,0,1])
        C = np.array([c1,c2],ndmin=2).T
        result = np.dot(A,x)

        self.assertEqual(True, np.array_equal(C,result))

    def test_zero_vector(self):
        a1 = Polynomial([0,1])
        a2 = Polynomial([0,0,1])
        a3 = Polynomial([1])
        a4 = Polynomial([0,1])

        b1 = Polynomial([0])
        b2 = Polynomial([0,0])

        A = np.array([[a1,a2],[a3,a4]], dtype=object)
        x = np.array([b1,b2],ndmin=2).T

        c1 = Polynomial([0,0,0,0])
        c2 = Polynomial([0,0,0])
        C = np.array([c1,c2],ndmin=2).T
        result = np.dot(A,x)

        self.assertEqual(True, np.array_equal(C,result))