import math
import numpy as np

def getN(epsilon,eta):
    N = np.array(np.zeros(4))
    N[0] = 0.25 * (1 - epsilon) * (1 - eta)
    N[1] = 0.25 * (1 + epsilon) * (1 - eta)
    N[2] = 0.25 * (1 + epsilon) * (1 + eta)
    N[3] = 0.25 * (1 - epsilon) * (1 + eta)

    return N


def getCordX(X,Y,px,py):
    N = getN(2* px-1,2* py-1)
    return X @ N

def getCordY(X,Y,px,py):
    N = getN(2*px - 1, 2* py -1)
    return Y @ N