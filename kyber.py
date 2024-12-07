from polynomial_helper import Polynomial
import numpy as np
from random import randrange
from decimal import Decimal, ROUND_HALF_UP

class KyberSystem():
    n: int
    k: int
    q: int

    def __init__(self, n, k, q, eta_1, eta_2):
        self.n = n
        self.k = k
        self.q = q
        self.eta_1 = eta_1
        self.eta_2 = eta_2

    def generate_key_pair(self):
        """

        :return: (A,t), s where (A,t) is the public key and s is the private key
        """
        s_polynomials = []
        for i in range(self.k):
            s_polynomials.append(Polynomial(self.generate_polynomial(self.n, self.q, self.eta_1)))
        s = np.array(s_polynomials, ndmin=2).T
        A = []
        for i in range(self.k):  #rows
            A_row = []
            for j in range(self.k):  # columns
                coeffs = []
                for coeff in range(self.n):
                    coeffs.append(randrange(self.q))
                A_row.append(Polynomial(coeffs))
            A.append(A_row)
        A = np.array(A, ndmin=2)
        e_polynomials = []
        for i in range(self.k):
            e_polynomials.append(Polynomial(self.generate_polynomial(self.n, self.q, self.eta_1)))
        e = np.array(e_polynomials, ndmin=2).T
        t = np.dot(A,s) + e
        return (A, t), s

    def encrypt(self, A, t, m):
        #r
        r_polynomials = []
        for i in range(self.k):
            r_polynomials.append(Polynomial(self.generate_polynomial(self.n, self.q, self.eta_1)))
        r = np.array(r_polynomials, ndmin=2).T
        #e1
        e1_polynomials = []
        for i in range(self.k):
            e1_polynomials.append(Polynomial(self.generate_polynomial(self.n, self.q, self.eta_2)))
        e1 = np.array(r_polynomials, ndmin=2).T
        #e2
        e2 = Polynomial(self.generate_polynomial(self.n, self.q, self.eta_2))

        m_polynomial = self.encode_message_to_polynomial(m)

        u = np.dot(A.T, r)+e1
        v = np.dot(t.T,r)+e2+m_polynomial
        return u, v

    def decrypt(self, s, u, v):
        message_noisy = v-(np.dot(s.T,u))
        message_noisy = message_noisy[0][0]
        message_noisy = message_noisy.coefficients.tolist()
        denoised = self.denoise(message_noisy)
        denoised_string = "".join([str(i) for i in denoised])
        converted_to_decimal = int(denoised_string, 2)
        return converted_to_decimal

    def denoise(self, coefficient_list):
        q = self.q
        q_half = self.get_nearest_integer_to_half(self.q)
        denoised_coefficient_list = []
        for coeff in coefficient_list:
            dist_q = abs(coeff-q)
            dist_0 = coeff
            dist_q_half = abs(coeff-q_half)
            if min(dist_q, dist_q_half, dist_0)==dist_q:
                denoised_coefficient_list.append(0)
            elif min(dist_q, dist_q_half, dist_0)==dist_0:
                denoised_coefficient_list.append(0)
            elif min(dist_q, dist_q_half, dist_0)==dist_q_half:
                denoised_coefficient_list.append(q_half)
        for i in range(len(denoised_coefficient_list)):
            denoised_coefficient_list[i] = int(denoised_coefficient_list[i]/q_half)
        return denoised_coefficient_list

    def get_nearest_integer_to_half(self, q):
        return int(Decimal(q/2).quantize(Decimal('1'), rounding=ROUND_HALF_UP))

    def encode_message_to_polynomial(self, m):
        m_binary_string = self.decimalToBinary(m)
        q_half = self.get_nearest_integer_to_half(self.q)
        coeffs = [int(i)*q_half for i in m_binary_string]
        poly = Polynomial(coeffs)
        return poly

    def decimalToBinary(self, n):
        return bin(n).replace("0b", "")
    def centered_binomial_eta(self, eta, n):
        """
        Generate coefficients sampled from the centered binomial distribution B_eta.

        Args:
            eta (int): Parameter defining the distribution, usually 2 or 3 in Kyber.
            n (int): Degree of the polynomial (number of coefficients), usually 256 in Kyber.

        Returns:
            np.ndarray: A numpy array of length n with coefficients sampled from B_eta.
        """
        # Generate 2 * eta * n random bits (0 or 1)
        random_bits = np.random.randint(0, 2, size=(2 * eta, n))

        # Split into two halves: a and b
        a = random_bits[:eta, :]  # First eta rows
        b = random_bits[eta:, :]  # Last eta rows

        # Compute the coefficients as the difference of sums
        coefficients = np.sum(a, axis=0) - np.sum(b, axis=0)

        return coefficients

    def generate_polynomial(self, n, q, eta):
        """
        Generate a polynomial with coefficients sampled from B_eta in R_q.

        Args:
            n (int): Degree of the polynomial (number of coefficients), usually 256.
            q (int): Modulus, usually 3329 in Kyber.
            eta (int): Parameter defining the binomial distribution, e.g., 2 or 3.

        Returns:
            np.ndarray: A polynomial of degree n with coefficients in Z_q.

        # Example usage
        n = 256  # Polynomial degree
        q = 3329  # Modulus
        eta = 2  # Parameter for binomial distribution

        poly = generate_polynomial(n, q, eta)
        print(poly)
        """
        # Generate coefficients from B_eta
        coefficients = self.centered_binomial_eta(eta, n)

        # Reduce coefficients modulo q
        polynomial = coefficients % q

        return polynomial


ks = KyberSystem(2, 4, 17, 3, 2)
(A,t), s = ks.generate_key_pair()
print(A)
print(t)
print(s)
print(ks.get_nearest_integer_to_half(3329))

u, v = ks.encrypt(A,t,24)
print(ks.encrypt(A,t,24))

decrypted = ks.decrypt(s,u,v)
print(decrypted)