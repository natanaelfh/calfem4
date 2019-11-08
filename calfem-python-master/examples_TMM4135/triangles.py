# -*- coding: utf-8 -*-
"""
Created on Sun Oct 21 08:15:51 2018

@author: bjohau
"""
import numpy as np


def plante(ex, ey, ep, D, eq=None):
    Dshape = D.shape
    if Dshape[0] != 3:
        raise NameError('Wrong constitutive dimension in plante')

    if ep[0] == 1:
        return tri3e(ex, ey, D, ep[1], eq)
    else:
        Dinv = np.inv(D)
        return tri3e(ex, ey, Dinv, ep[1], eq)


def tri3e(ex, ey, D, th, eq=None):
    """
    Compute the stiffness matrix for a two dimensional beam element.
    
    :param list ex: element x coordinates [x1, x2, x3]
    :param list ey: element y coordinates [y1, y2, y3]
    :param list D : 2D constitutive matrix
    :param list th: element thickness
    :param list eq: distributed loads, local directions [bx, by]
    :return mat Ke: element stiffness matrix [6 x 6]
    :return mat fe: consistent load vector [6 x 1] (if eq!=None)
    """

    tmp = np.array([[1, ex[0], ey[0]],
                    [1, ex[1], ey[1]],
                    [1, ex[2], ey[2]]])

    A2 = np.linalg.det(tmp)  # Double of triangle area
    A = A2 / 2.0

    cyclic_ijk = [0, 1, 2]  # Cyclic permutation of the nodes i,j,k

    zi_px, zi_py = zeta_partials_x_and_y(ex, ey)

    B = np.array([
        [zi_px[0], 0, zi_px[1], 0, zi_px[2], 0],
        [0, zi_py[0], 0, zi_py[1], 0, zi_py[2]],
        [zi_py[0], zi_px[0], zi_py[1], zi_px[1], zi_py[2], zi_px[2]]])

    Ke = B.T @ D @ B * A * th

    if eq is None:
        return Ke
    else:
        fx = A * th * eq[0] / 3.0
        fy = A * th * eq[1] / 3.0
        fe = np.array([[fx], [fy], [fx], [fy], [fx], [fy]])
        return Ke, fe


def zeta_partials_x_and_y(ex, ey):
    """
    Compute partials of area coordinates with respect to x and y.
    
    :param list ex: element x coordinates [x1, x2, x3]
    :param list ey: element y coordinates [y1, y2, y3]
    """

    tmp = np.array([[1, ex[0], ey[0]],
                    [1, ex[1], ey[1]],
                    [1, ex[2], ey[2]]])

    A2 = np.linalg.det(tmp)  # Double of triangle area

    cyclic_ijk = [0, 1, 2]  # Cyclic permutation of the nodes i,j,k

    zeta_px = np.zeros(3)  # Partial derivative with respect to x
    zeta_py = np.zeros(3)  # Partial derivative with respect to y

    for i in cyclic_ijk:
        j = (i + 1) % 3
        k = (i + 2) % 3

        zeta_px[i] = (ey[j] - ey[k]) / A2
        zeta_py[i] = (ex[k] - ex[j]) / A2

    return zeta_px, zeta_py


# Functions for 6 node triangle


def tri6_area(ex, ey):
    tmp = np.array([[1, ex[0], ey[0]],
                     [1, ex[1], ey[1]],
                     [1, ex[2], ey[2]]])

    A = np.linalg.det(tmp) / 2

    return A


def tri6_shape_functions(zeta):
    cyclic_ijk = [0, 1, 2, 3, 4, 5]  # Cyclic permutation of the nodes i,j,k
    L = getL(ex)
    N6 = np.zeros(6)
    N6[0] = 6

    # TODO: fill out missing parts (or reformulate completely)

    return N6

def getL(ex,ey,x,y):
    cyclic = [0,1,2]
    constants = np.array(np.zeros((3,3)))
    A = tri6_area(ex,ey)
    for i in range(3):
        j = cyclic[(i+1)%3]
        k = cyclic[(i+2)%3]
        constants[i,0] = ex[j] * ey[k] - ex[k] * ey[j]
        constants[i,1] = ey[j] - ey[k]
        constants[i,2] = ex[k] - ex[j]

    L = np.array(np.zeros((3,1)))
    for i in range(3):
        a = constants[i,0]
        b = constants[i,1]
        c = constants[i,2]

        L[i] = (a + b * x + c * y)/(2*A)

    return L


def tri6_shape_function_partials_x_and_y(ex, ey, x, y):
    zeta_px, zeta_py = zeta_partials_x_and_y(ex, ey)

    N6_px = np.zeros(6)
    N6_py = np.zeros(6)

    cyclic_ijk = [0, 1, 2, 3, 4, 5]  # Cyclic permutation of the nodes i,j,k

    # TODO: fill out missing parts (or reformulate completely)

    return N6_px, N6_py


def tri6_Bmatrix(zeta, ex, ey):
    nx, ny = tri6_shape_function_partials_x_and_y(zeta, ex, ey)

    Bmatrix = np.matrix(np.zeros((3, 12)))

    # TODO: fill out missing parts (or reformulate completely)

    return Bmatrix


def tri6_Kmatrix(ex, ey, D, th, eq=None):
    zetaInt = np.array([[0.5, 0.5, 0.0],
                        [0.0, 0.5, 0.5],
                        [0.5, 0.0, 0.5]])

    wInt = np.array([1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0])

    A = tri6_area(ex, ey)

    Ke = np.array(np.zeros((12, 12)))

    # TODO: fill out missing parts (or reformulate completely)

    if eq is None:
        return Ke
    else:
        fe = np.array(np.zeros((12, 1)))

        # TODO: fill out missing parts (or reformulate completely)

        return Ke, fe


def tri6e(ex, ey, D, th, eq=None):
    return tri6_Kmatrix(ex, ey, D, th, eq)
