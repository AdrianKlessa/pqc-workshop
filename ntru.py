import numpy as np
from polynomial_helper import Polynomial
def get_ternary_polynomial(no_positive,no_negative, degree):
    """
    :param no_positive: Number of 1s in the vector
    :param no_negative: Number of -1s in the vector
    :param degree: Maximum degree of the vector
    :return:
    """
    vec = np.zeros(degree+1)
    indices_positive = np.random.choice(degree+1, no_positive, replace=False)
    vec[indices_positive] = 1
    possible_negative_indices = np.arange(degree+1)
    possible_negative_indices = np.delete(possible_negative_indices, indices_positive)
    indices_negative = np.random.choice(possible_negative_indices, no_negative, replace=False)
    vec[indices_negative] = -1
    return Polynomial(vec)

def center_coefficients(coefficients, q):
    """
    Center the coefficients of a polynomial with respect to q, such that for each coefficient a
    q/2 <= a < q
    :param coefficients: Coefficients to center
    :param q: Center value
    :return: Coefficients centered to range [-q/2,q/2)
    """
    centered = []
    for c in coefficients:
        if c>=q/2:
            centered.append((-q/2)+(c%(q/2)))
        elif c<-q/2:
            centered.append((c%(q/2)))
        else:
            centered.append(c)
    return centered
