from polynomial_helper import Polynomial
import numpy as np
from random import randrange
from decimal import Decimal, ROUND_HALF_UP


class KyberSystem():
    n: int
    k: int
    q: int
    modulo_polynomial: Polynomial

    def __init__(self, n=256, k=4, q=3329, eta_1=2, eta_2=2):
        """
        Create an instance of the Kyber key encapsulation mechanism
        :param n: Maximum degree of polynomials used
        :param k: Number of polynomials in a vector. Arrays have dimension k x k
        :param q: Modulus for coefficients of polynomials
        :param eta_1: Eta1 parameter for generation of error vectors
        :param eta_2: Eta2 parameter for generation of error vectors
        """
        self.n = n
        self.k = k
        self.q = q
        self.eta_1 = eta_1
        self.eta_2 = eta_2
        modulo_polynomial_coeffs = [1] + ([0] * q)
        modulo_polynomial_coeffs[-1] = 1
        self.modulo_polynomial = Polynomial(modulo_polynomial_coeffs)

    def reduce_array_modulo(self, array):
        """
        Reduce every element of the input array (coefficients reduced mod q, polynomail reduced mod (X^n)+1)
        :param array: Array to be reduced
        :return: Input array, with each element being a reduced polynomial
        """
        for (i, j), value in np.ndenumerate(array):
            array[i][j] = self.reduce_element_modulo(array[i][j])

    def reduce_element_modulo(self, p: Polynomial):
        """

        :param p: A polynomial
        :return: The polynomial p with coefficients reduced modulo q, itself reduced modulo the polynomial (X^n)+1
        """
        a = p.reduced_modulo_scalar(self.q)
        _, a = a.divide_by(self.modulo_polynomial, self.q)
        return a

    def generate_key_pair(self):
        """
        Generates the Kyber key pair
        :return: (A,t), s where (A,t) is the public key and s is the private key
        """
        # Generate secret key
        s_polynomials = []
        for i in range(self.k):
            s_polynomials.append(Polynomial(self.generate_polynomial(self.n, self.q, self.eta_1)))
        s = np.array(s_polynomials, ndmin=2).T

        # Generate matrix A for the public key
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

        # Generate e
        e_polynomials = []
        for i in range(self.k):
            e_polynomials.append(Polynomial(self.generate_polynomial(self.n, self.q, self.eta_1)))
        e = np.array(e_polynomials, ndmin=2).T

        # Calculate t for the public key
        t = np.dot(A, s) + e
        self.reduce_array_modulo(t)
        return (A, t), s

    def encrypt(self, A, t, m: int):
        """
        Encrypt message m using the public key A, t
        :param A: Matrix part of public key
        :param t: Vector part of public key
        :param m: Message to be encrypted (as int)
        :return: Encrypted message (u,v)
        """
        #r
        r_polynomials = []
        for i in range(self.k):
            r_polynomials.append(Polynomial(self.generate_polynomial(self.n, self.q, self.eta_1)))
        r = np.array(r_polynomials, ndmin=2).T

        #e1
        e1_polynomials = []
        for i in range(self.k):
            e1_polynomials.append(Polynomial(self.generate_polynomial(self.n, self.q, self.eta_2)))
        e1 = np.array(e1_polynomials, ndmin=2).T

        #e2
        e2 = Polynomial(self.generate_polynomial(self.n, self.q, self.eta_2))

        m_polynomial = self.encode_message_to_polynomial(m)
        u = np.dot(A.T, r) + e1
        self.reduce_array_modulo(u)
        v = np.dot(t.T, r) + e2 + m_polynomial
        self.reduce_array_modulo(v)
        return u, v

    def decrypt(self, s, u, v) -> int:
        """
        Decrypt message u,v using the public key s
        :param s: Secret key
        :param u: u part of the public key
        :param v: v part of the public key
        :return: Decrypted message (int)
        """
        message_noisy = v - (np.dot(s.T, u))
        message_noisy = message_noisy[0][0]
        message_noisy = self.reduce_element_modulo(message_noisy)
        message_noisy = message_noisy.coefficients.tolist()

        denoised = self.denoise(message_noisy)
        denoised_string = "".join([str(i) for i in denoised])
        # Reverse the reversion from polynomial encoding using [::-1]
        converted_to_decimal = int(denoised_string.rstrip("0")[::-1], 2)
        return converted_to_decimal

    def denoise(self, coefficient_list):
        """
        Returns 0 for every coefficient in coefficient_list closer to q or 0, 1 for every coefficient closer to q/2
        :param coefficient_list: Coefficents of the noisy message
        :return: Denoised coefficients
        """
        q = self.q
        q_half = self.get_nearest_integer_to_half(self.q)
        denoised_coefficient_list = []
        for coeff in coefficient_list:
            dist_q = abs(coeff - q)
            dist_0 = coeff
            dist_q_half = abs(coeff - q_half)
            if min(dist_q, dist_q_half, dist_0) == dist_q:
                denoised_coefficient_list.append(0)
            elif min(dist_q, dist_q_half, dist_0) == dist_0:
                denoised_coefficient_list.append(0)
            else:
                denoised_coefficient_list.append(q_half)
        for i in range(len(denoised_coefficient_list)):
            denoised_coefficient_list[i] = int(denoised_coefficient_list[i] / q_half)
        return denoised_coefficient_list

    def get_nearest_integer_to_half(self, q: int) -> int:
        """
        Custom function for finding the nearest integer to q/2, rounding ties up as defined in the Kyber specification
        :param q:
        :return: Number closest to q/2, with ties rounded up
        """
        return int(Decimal(q / 2).quantize(Decimal('1'), rounding=ROUND_HALF_UP))

    def encode_message_to_polynomial(self, m: int) -> Polynomial:
        """
        Encode a given message (int) into a polynomial with binary coefficients
        :param m: Message to encode
        :return: Encoded message (Polynomial)
        """
        m_binary_string = self.decimalToBinary(m)
        q_half = self.get_nearest_integer_to_half(self.q)

        # Reverse to prevent truncation of messages ending with 0
        coeffs = list(reversed([int(i) * q_half for i in m_binary_string]))
        poly = Polynomial(coeffs)
        return poly

    def decimalToBinary(self, n):
        """
        Custom function for converting a decimal number to a binary string without the 0b prefix
        :param n:
        :return:
        """
        return bin(n).replace("0b", "")

    def get_coeff_from_eta(self, eta):
        """
        Get a random coefficient from the Centered Binomial Distribution defined by parameter eta
        :param eta: Value of eta, as defined in the Kyber specification
        :return: A single number from the eta Centered Binomial Distribution
        """
        random_bits_1 = np.random.randint(0, 2, size=eta)
        random_bits_2 = np.random.randint(0, 2, size=eta)
        return abs(random_bits_1.sum() - random_bits_2.sum())

    def generate_polynomial(self, n, q, eta):
        """
        Generate a polynomial with n coefficients from the Centered Binomial Distribution defined by parameter eta
        :param n: Number of coefficients to generate
        :param q: Modulus for coefficients
        :param eta: Parameter for the Centered Binomial Distribution
        :return:
        """
        return [self.get_coeff_from_eta(eta) % q for _ in range(n)]
