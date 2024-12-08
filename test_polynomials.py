import unittest
from polynomial_helper import Polynomial, multiplicative_inverse


class TestPolynomials(unittest.TestCase):
    def test_equality(self):
        p1 = Polynomial([0, 0, 1])
        p2 = Polynomial([0, 0, 1])

        self.assertEqual(p1, p2)

    def test_inequality(self):
        p1 = Polynomial([0, 0, 1])
        p2 = Polynomial([1, 0, 0, 0])
        self.assertNotEqual(p1, p2)

    def test_coefficient_cutoff_1(self):
        p1 = Polynomial([1, 1, 0])
        p2 = Polynomial([1, 1])

        self.assertEqual(p1, p2)

    def test_coefficient_cutoff_2(self):
        p1 = Polynomial([0,0,0])
        p2 = Polynomial([0])

        self.assertEqual(p1, p2)
        self.assertEqual(len(p1.coefficients),1)

    def test_addition_zero(self):
        p1 = Polynomial([0, 1])
        p2 = Polynomial([0])

        a1 = p1.add_mod(p2, 3)
        self.assertEqual(a1, p1)

    def test_multiplication_one(self):
        p1 = Polynomial([1, 2, 3, 4])
        p2 = Polynomial([1])

        a1 = p1.multiply_mod(p2, 7, None)
        self.assertEqual(a1, p1)

    def test_division_zero(self):
        p1 = Polynomial([1, 2, 3, 4])
        p2 = Polynomial([0])
        method_to_test = p1.divide_by
        self.assertRaises(ValueError, method_to_test, p2, None)

    def test_add(self):
        p1 = Polynomial([0, 0, 1])
        p2 = Polynomial([1, 0, 0, 0])
        a1 = p1.add_mod(p2, None)
        self.assertEqual(a1, Polynomial([1, 0, 1, 0]))
        self.assertEqual(p1+p2, Polynomial([1, 0, 1, 0]))

        p3 = Polynomial([0, 0, 1])
        p4 = Polynomial([0, 0, 4])
        a2 = p3.add_mod(p4, None)
        self.assertEqual(a2, Polynomial([0, 0, 5]))

    def test_add_modulo(self):
        p1 = Polynomial([0, 0, 1])
        p2 = Polynomial([1, 0, 0, 1])
        a1 = p1.add_mod(p2, 2)
        self.assertEqual(a1, Polynomial([1, 0, 1, 1]))

        p3 = Polynomial([0, 0, 1])
        p4 = Polynomial([0, 0, 4])
        a2 = p3.add_mod(p4, 2)
        self.assertEqual(a2, Polynomial([0, 0, 1]))

    def test_sub(self):
        p1 = Polynomial([0, 0, 1])
        p2 = Polynomial([1, 0, 0, 1])
        a1 = p1.substract_mod(p2, None)
        self.assertEqual(a1, Polynomial([-1, 0, 1, -1]))
        self.assertEqual(p1-p2, Polynomial([-1, 0, 1, -1]))

    def test_sub_modulo(self):
        p1 = Polynomial([0, 0, 1])
        p2 = Polynomial([3, 0, 0, 3])
        a1 = p1.substract_mod(p2, 2)
        self.assertEqual(a1, Polynomial([1, 0, 1, 1]))

    def test_reduce_modulo(self):
        p1 = Polynomial([0, 0, 5]).reduced_modulo_scalar(7)
        p2 = Polynomial([14, 0, 8]).reduced_modulo_scalar(7)
        self.assertEqual(p1, Polynomial([0, 0, 5]))
        self.assertEqual(p2, Polynomial([0, 0, 1]))

    def test_mul(self):
        p1 = Polynomial([0, 0, 2])
        p2 = Polynomial([1, -5, 3])
        a1 = p1.multiply_mod(p2, None, None)
        self.assertEqual(a1, Polynomial([0, 0, 2, -10, 6]))
        self.assertEqual(p1*p2, Polynomial([0, 0, 2, -10, 6]))

    def test_mul_modulo(self):
        p1 = Polynomial([2, 3, 0, 0, 10])
        p2 = Polynomial([13, 2, 0, 1])
        a1 = p1.multiply_mod(p2, 17, None)
        self.assertEqual(a1, Polynomial([9, 9, 6, 2, 14, 3, 0, 10]))

    def test_add_zero(self):
        p1 = Polynomial([0, 0, 1])
        p2 = Polynomial([0])
        a1 = p1.add_mod(p2, None)
        self.assertEqual(a1, p1)

    def test_subtract_zero(self):
        p1 = Polynomial([0, 0, 1])
        p2 = Polynomial([0])
        a1 = p1.substract_mod(p2, None)
        self.assertEqual(a1, p1)

    def test_multiply_zero(self):
        p1 = Polynomial([0, 0, 1])
        p2 = Polynomial([0])
        a1 = p1.multiply_mod(p2, None, None)
        self.assertEqual(a1, p2)

    def test_mul_modulo_polynomial(self):
        a = Polynomial([1, 1, 1])  # 1+x+x^2
        b = Polynomial([0, 1, 0, 1])  # x+x^3
        mod = Polynomial([1, 0, 0, 0, 1])  # 1+x^4
        c = a.multiply_mod(b, 2, mod)  # Should be x^2+1
        self.assertEqual(c, Polynomial([1, 0, 1]))

    def test_div(self):
        p1 = Polynomial([-4, 0, -2, 1])
        p2 = Polynomial([-3, 1])

        q, r = p1.divide_by(p2, None)
        self.assertEqual(q, Polynomial([3, 1, 1]))
        self.assertEqual(r, Polynomial([5]))

    def test_div_modulo(self):
        p1 = Polynomial([9, 9, 6, 2, 14, 3, 0, 10])
        p2 = Polynomial([13, 2, 0, 1])
        a1, _ = p1.divide_by(p2, 17)
        self.assertEqual(a1, Polynomial([2, 3, 0, 0, 10]))

        p3 = Polynomial([-1, 0, -1, 1])
        p4 = Polynomial([1, 1])

        a1, _ = p3.divide_by(p4, 3)
        self.assertEqual(a1, Polynomial([2, 1, 1]))

        p5 = Polynomial([-1,0,-1,1])
        p6 = Polynomial([-5,1])
        a1, _ = p5.divide_by(p6, 11)
        self.assertEqual(a1, Polynomial([9,4,1]))

    def test_division_self(self):
        p1 = Polynomial([9, 9, 6, 2, 14, 3, 0, 10])
        p2 = Polynomial([9, 9, 6, 2, 14, 3, 0, 10])
        a1, _ = p1.divide_by(p2, 17)
        self.assertEqual(a1, Polynomial([1]))


    def test_div_modulo_remainder(self):
        p1 = Polynomial([7, 10, 5, 2])
        p2 = Polynomial([4, 0, 1])
        r = Polynomial([4, 2])
        q = Polynomial([5, 2])
        a1, a2 = p1.divide_by(p2, 17)

        self.assertEqual(a1, q)
        self.assertEqual(a2, r)

    def test_inverse_1(self):
        a = Polynomial([1, 0, 1])  # x^2 + 1
        p = Polynomial([1, 2, 0, 1])  # x^3+2x+1
        a1 = a.get_inverse(p, 3)

        self.assertEqual(a1, Polynomial([2, 1, 2]))

    def test_inverse_2(self):
        a = Polynomial([1, 1, 0, 0, 1, 0, 1])  # x^6 +x^4+x+1
        p = Polynomial([1, 1, 0, 1, 1, 0, 0, 0, 1])  # x^8+x^4+x^3+x+1
        a1 = a.get_inverse(p, 2)

        self.assertEqual(a1, Polynomial([0, 1, 0, 1, 0, 0, 1, 1]))

    def test_inverse_3(self):
        a = Polynomial([0, 1])
        p = Polynomial([1, 4, 1])

        a1 = a.get_inverse(p, 7)
        self.assertEqual(a1, Polynomial([3, 6]))

    def test_inverse_4(self):
        a = Polynomial([0, 1, 0, 1])
        p = Polynomial([1, 1, 0, 0, 1])
        a1 = a.get_inverse(p, 2)

        self.assertEqual(a1, Polynomial([0, 0, 1, 1]))

    def test_inverse_5(self):
        a = Polynomial([1,1,0,0,1])
        p = Polynomial([-1,0,0,0,0,1])
        a1 = a.get_inverse(p, 2)
        self.assertEqual(a1, Polynomial([1,0,1,1]))

    def test_inverse_ntru(self):
        N=11
        f = Polynomial([-1, 0, 1, -1,0,0,0,0,1,0,1])
        q = 73

        N_polynomial = Polynomial([-1,0,0,0,0,0,0,0,0,0,0,1])
        f_q = f.get_inverse(N_polynomial, q)
        print(f)
        print(N_polynomial)
        self.assertEqual(f_q, Polynomial([28,-32,12,-30,-33,36,-10,33,15,33,22]).reduced_modulo_scalar(q))

    def test_noninvertible_polynomial(self):
        a = Polynomial([0, 1])
        p = Polynomial([0, 1, 1])

        method_to_test = a.get_inverse
        self.assertRaises(ValueError, method_to_test, p, 2)

    def test_inverse_scalar(self):
        a1 = multiplicative_inverse(3, 11)
        a2 = multiplicative_inverse(2, 11)
        a3 = multiplicative_inverse(78, 155)

        self.assertEqual(a1, 4)
        self.assertEqual(a2, 6)
        self.assertEqual(a3, 2)
