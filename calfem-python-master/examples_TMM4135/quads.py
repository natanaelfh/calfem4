# -*- coding: utf-8 -*-
"""
Created on Sun Oct 21 08:15:51 2018

@author: bjohau
"""
import numpy as np


def quad4e(ex, ey, D, th, eq=None):
    """
    Compute the stiffness matrix for a four node membrane element.

    :param list ex: element x coordinates [x1, x2, x3]
    :param list ey: element y coordinates [y1, y2, y3]
    :param list D : 2D constitutive matrix
    :param list th: element thickness
    :param list eq: distributed loads, local directions [bx, by]
    :return mat Ke: element stiffness matrix [6 x 6]
    :return mat fe: consistent load vector [6 x 1] (if eq!=None)
    """

    ex1 = ex[0:3]
    ex2 = ex[1:4]
    ey1 = ey[0:3]
    ey2 = ey[1:4]

    Ke = np.array(np.zeros((8, 8)))
    fe = np.array(np.zeros((8, 1)))

    """"
    MY CODE:
    
    
    
    
    """
    '''
    yDn = (ey[0] - ey[3]) * 0.5
    yDe = (ey[0] - ey[1]) * 0.5
    xDn = (ex[0] - ex[3]) * 0.5
    xDe = (ex[0] - ex[1]) * 0.5
    NDe = np.array([0.5, -0.5, 0, 0])
    NDn = np.array([0.5, 0, 0, -0.5])

    A = xDe * yDn - xDn * yDe

    B = np.array(np.zeros((2, 4)))
    for i in range(4):
        B[0, i] = (yDn * NDe[i] - yDe * NDn[i]) / A
        B[1, i] = (-xDn * NDe[i] + xDe * NDn[i]) / A

    tmp1 = np.array([[1, ex[0], ey[0]],
                     [1, ex[1], ey[1]],
                     [1, ex[2], ey[2]]])

    tmp2 = np.array([[1, ex[1], ey[1]],
                     [1, ex[2], ey[2]],
                     [1, ex[3], ey[3]]])

    A1 = np.linalg.det(tmp1)
    A2 = np.linalg.det(tmp2)

    A = A1 / 2 + A2 / 2
    BB = np.zeros((3, 8))

    for i in range(3):
        if i == 0:
            for j in range(4):
                BB[i, j] = B[0, j]
        if i == 1:
            for j in range(4):
                BB[i, j + 4] = B[1, j]
        if i == 2:
            for j in range(4):
                BB[i, j] = B[1, j]
                BB[i, j + 4] = B[0, j]

    Ke = BB.T @ D @ BB * A * th
    '''
    '''
    Nytt forsøk under
    '''
    tmp1 = np.array([[1, ex[0], ey[0]],
                     [1, ex[1], ey[1]],
                     [1, ex[2], ey[2]]])

    tmp2 = np.array([[1, ex[2], ey[2]],
                     [1, ex[3], ey[3]],
                     [1, ex[0], ey[0]]])

    A1 = np.linalg.det(tmp1)
    A2 = np.linalg.det(tmp2)

    A = A1 / 2 + A2 / 2

    B = np.array(np.zeros((3, 8)))
    Ke = np.array(np.zeros((8, 8)))
    for i in range(-1,1,1):
        for j in range(-1,1,1):
            B = getB(ex, ey, i*0.7745966, j*0.7745966)
            J  = getJacobi(ex,ey , i*0.7745966, j*0.7745966)
            Jdet = np.linalg.det(J)
            if i == 0:
                H1 = 0.888888
            else:
                H1 = 0.555555
            if j == 0:
                H2 = 0.888888
            else:
                H2 = 0.555555
            Ke = Ke + (B.T @ D @ B) * Jdet * th * H1 * H2 * 1

    # TODO: fill out missing parts (or reformulate completely)
    # Bruk numerisk integrasjon
    cyclic = [0, 1, 2, 3]

    if eq is None:
        return Ke
    else:
        fe = np.array(np.zeros((8,1)))
        for i in range(4):
            j = cyclic[(i + 1) % 4]
            k = cyclic[(i + 3) % 4]
            fe[2*i] = eq[0] * th * 2
            fe[2*i + 1] = eq[1] * th * 2

        return Ke, fe

def getN(ex, ey, epsilon, eta):
    N = np.array(np.zeros(1, 4))
    N[0] = 0.25 * (1+epsilon)*(1+eta)
    N[1] = 0.25 * (1-epsilon)*(1+eta)
    N[2] = 0.25 * (1-epsilon)*(1-eta)
    N[3] = 0.25 * (1+epsilon)*(1-eta)



def getN_eta(ex, ey, epsilon):
    N_eta = np.array(np.zeros(4))
    N_eta[0] = 0.25 * (1 + epsilon)
    N_eta[1] = 0.25 * (1 - epsilon)
    N_eta[2] = -0.25 * (1 - epsilon)
    N_eta[3] = -0.25 * (1 + epsilon)

    return N_eta


def getN_epsilon(ex, ey, eta):
    N_epsilon = np.array(np.zeros(4))
    N_epsilon[0] = 0.25 * (1 + eta)
    N_epsilon[1] = -0.25 * (1 + eta)
    N_epsilon[2] = -0.25 * (1 - eta)
    N_epsilon[3] = 0.25 * (1 - eta)

    return N_epsilon


def gety_eta(ex, ey, epsilon):
    N_eta = getN_eta(ex, ey, epsilon)
    y_eta = 0.0
    for i in range(len(ey)):
        y_eta += N_eta[i] * ey[i]
    return y_eta


def gety_epsilon(ex, ey, eta):
    N_epsilon = getN_epsilon(ex, ey, eta)
    y_epsilon = 0.0
    for i in range(len(ey)):
        y_epsilon += N_epsilon[i] * ey[i]
    return y_epsilon


def getx_eta(ex, ey, epsilon):
    N_eta = getN_eta(ex, ey, epsilon)
    x_eta = 0.0
    for i in range(len(ex)):
        x_eta += N_eta[i] * ex[i]
    return x_eta


def getx_epsilon(ex, ey, eta):
    N_epsilon = getN_epsilon(ex, ey, eta)
    x_epsilon = 0.0
    for i in range(len(ex)):
        x_epsilon += N_epsilon[i] * ex[i]
    return x_epsilon


def getN_x(ex, ey, epsilon, eta):
    N_x = np.array(np.zeros(4))
    y_eta = gety_eta(ex, ey, epsilon)
    y_epsilon = gety_epsilon(ex, ey, eta)
    x_eta = getx_eta(ex, ey, epsilon)
    x_epsilon = getx_epsilon(ex, ey, eta)

    N_epsilon = getN_epsilon(ex, ey, eta)
    N_eta = getN_eta(ex, ey, epsilon)

    for i in range(len(ex)):
        N_x[i] = y_eta * N_epsilon[i] - y_epsilon * N_eta[i]

    return N_x / (x_epsilon * y_eta - x_eta * y_epsilon)


def getN_y(ex, ey, epsilon, eta):
    N_y = np.array(np.zeros(4))
    y_eta = gety_eta(ex, ey, epsilon)
    y_epsilon = gety_epsilon(ex, ey, eta)
    x_eta = getx_eta(ex, ey, epsilon)
    x_epsilon = getx_epsilon(ex, ey, eta)

    N_epsilon = getN_epsilon(ex, ey, eta)
    N_eta = getN_eta(ex, ey, epsilon)

    for i in range(len(ey)):
        N_y[i] = x_epsilon * N_eta[i] - x_eta * N_epsilon[i]

    return N_y / (x_epsilon * y_eta - x_eta * y_epsilon)


def getB(ex,ey,epsilon,eta):
    B = np.array(np.zeros((3, 2 * len(ex))))
    N_x = getN_x(ex, ey, epsilon, eta)
    N_y = getN_y(ex, ey, epsilon, eta)

    for i in range(3):
        if i == 0:
            for j in range(4):
                B[i, 2*j] = N_x[j]
        if i == 1:
            for j in range(4):
                B[i, 2*j +1] = N_y[j]
        if i == 2:
            for j in range(4):
                B[i, 2*j] = N_y[j]
                B[i, 2*j +1] = N_x[j]
    return B


def getJacobi(ex,ey, epsilon, eta):
    x_epsilon = getx_epsilon(ex,ey,eta)
    x_eta = getx_eta(ex, ey, epsilon)
    y_epsilon = gety_epsilon(ex,ey, eta)
    y_eta = gety_eta(ex, ey, epsilon)

    jacobi = np.array([[x_epsilon, y_epsilon],[x_eta, y_eta]])

    return jacobi


def quad9e(ex, ey, D, th, eq=None):
    """
    Compute the stiffness matrix for a four node membrane element.

    :param list ex: element x coordinates [x1, x2, x3]
    :param list ey: element y coordinates [y1, y2, y3]
    :param list D : 2D constitutive matrix
    :param list th: element thickness
    :param list eq: distributed loads, local directions [bx, by]
    :return mat Ke: element stiffness matrix [6 x 6]
    :return mat fe: consistent load vector [6 x 1] (if eq!=None)
    """

    Ke = np.matrix(np.zeros((18, 18)))
    fe = np.matrix(np.zeros((18, 1)))

    # TODO: fill out missing parts (or reformulate completely)

    if eq is None:
        return Ke
    else:
        return Ke, fe

def funct(x,y):
    return x*y

def dIntegral(f,x_0,x_1,y_0,y_1):
    return f(x_0,y_0)
