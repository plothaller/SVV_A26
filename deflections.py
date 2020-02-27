import numpy as np
import Aerodynamic_Loading_Main_V3 as AV3

#The deflections.py file contains the functions that produce the different deflection plots

def macaulay(x, x_n, power):
  result = (x-x_n)
  if result >= 0:
    return result**power
  else:
    return 0

def internal_moment_y(la, F_I, F_1z, F_2z, F_3z, x1, x2, x3, xa, theta, P):
    # Temporary variables:
    xi1 = x2 - xa/2
    xi2 = x2 + xa/2
    x = np.linspace(0, la, 500)
    y = []

    for i in x:
        y.append(-F_1z*macaulay(i, x1, 1) - F_I * np.cos(theta) * macaulay(i, xi1, 1) - F_2z*macaulay(i, x2, 1) + P *np.cos(theta)*macaulay(i, xi2, 1) - F_3z*macaulay(i, x3, 1))

    return x, y


def internal_moment_z(la, F_I, F_1y, F_2y, F_3y, x1, x2, x3, xa, theta, P, zsc):
    # Temporary variables:
    xi1 = x2 - xa/2
    xi2 = x2 + xa/2
    x = np.linspace(0, la, 500)
    y = []


    for i in x:
        y.append(-F_1y*macaulay(i, x1, 1) - F_I*np.sin(theta)*macaulay(i, xi1, 1) - F_2y * macaulay(i, x2, 1) + P*np.sin(theta)*macaulay(i, xi2, 1) - F_3y*macaulay(i, x3, 1) + AV3.ThreeIntegral(i))

    return x, y

def torque(la, F_I, F_1y, F_2y, F_3y, x1, x2, x3, xa, theta, P, zsc, ha):
    # Temporary variables:
    xi1 = x2 - xa / 2
    xi2 = x2 + xa / 2
    x = np.linspace(0, la, 500)
    y = []


    for i in x:
        y.append((zsc*F_1y*macaulay(i, x1, 0) + (ha/2)*F_I*np.cos(theta)*macaulay(i, xi1, 0) + zsc*F_I*np.sin(theta)*macaulay(i, xi1, 0) + zsc* F_2y*macaulay(i, x2, 0) - zsc*P*np.sin(theta)*macaulay(i, xi2, 0) + zsc*F_3y*macaulay(i, x3, 0) +AV3.DoubleIntegralZSC(i, zsc)))

    return x, y

def v_deflection(la, F_I, F_1y, F_2y, F_3y, x1, x2, x3, xa, theta, P, c1, c2,  E, I_zz):
    #Printing without aeroload
    # Temporary variables:
    xi1 = x2 - xa/2
    xi2 = x2 + xa/2
    k = -1/(E*I_zz)
    x_v = np.linspace(0, la, 500)
    y_v = []


    for i in x_v:
        y_v.append(k*(-(F_1y/6)*macaulay(i, x1, 3) - (F_I/6)*np.sin(theta)*macaulay(i, xi1, 3) - (F_2y/6)*macaulay(i, x2, 3) + (P/6)*np.sin(theta)*macaulay(i, xi2, 3) - (F_3y/6)*macaulay(i, x3, 3) + 0) + c1*i + c2) #AV3.FiveIntegral(i)

    return x_v, y_v

def w_deflection(la, F_I, F_1z, F_2z, F_3z, x1, x2, x3, xa, theta, P, c3, c4, E, I_yy):
    # Printing without aeroload
    # Temporary variables:
    xi1 = x2 - xa/2
    xi2 = x2 + xa/2
    k = -1 / (E * I_yy)
    x_w = np.linspace(0, la, 500)
    y_w = []

    for i in x_w:
        y_w.append((k * (-(F_1z / 6) * macaulay(i, x1, 3) - (F_I / 6) * np.cos(theta) * macaulay(i, xi1, 3) - (F_2z / 6) * macaulay(i, x2, 3) + (P / 6) * np.cos(theta) * macaulay(i, xi2, 3) - (F_3z / 6) * macaulay(i, x3, 3)) + c3 * i + c4))

    return x_w, y_w

def shear_z(la, F_I, F_1z, F_2z, F_3z, x1, x2, x3, xa, theta, P):
    #Printing without aeroload
    # Temporary variables:
    xi1 = x2 - xa/2
    xi2 = x2 + xa/2
    x = np.linspace(0, la, 500)
    y = []

    for i in x:
        y.append(-F_1z*macaulay(i, x1, 0) - F_I*np.cos(theta)*macaulay(i, xi1, 0) - F_2z*macaulay(i, x2, 0) + P*np.cos(theta)*macaulay(i, xi2, 0) - F_3z*macaulay(i, x3, 0))

    return x, y

def shear_y(la, F_I, F_1y, F_2y, F_3y, x1, x2, x3, xa, theta, P):
    #Printing without aeroload
    # Temporary variables:
    xi1 = x2 - xa/2
    xi2 = x2 + xa/2
    x = np.linspace(0, la, 500)
    y = []
    DOUBLEINTEGRAL = 0

    for i in x:
        y.append(-F_1y*macaulay(i, x1, 0) - F_I*np.sin(theta)*macaulay(i, xi1, 0) - F_2y*macaulay(i, x2, 0) + P*np.sin(theta)*macaulay(i, xi2, 0) - F_3y*macaulay(i, x3, 0) + AV3.DoubleIntegral(i))

    return x, y


def Twist(la, F_I, F_1y, F_2y, F_3y, x1, x2, x3, xa, theta, P, c5, G, J, zsc, ha):
    xi1 = x2 - xa / 2  # points where the actuators act
    xi2 = x2 + xa / 2
    k = 1 / (G * J)
    x = np.linspace(0, la, 500)
    y = []


    for i in x:
        y.append((k*((F_1y*zsc)*macaulay(i, x1, 1) + (F_2y*zsc)*macaulay(i, x2, 1) + (F_I*(ha/2)*np.cos(theta))*macaulay(i, xi1, 1) + (F_I*zsc*np.sin(theta))*macaulay(i, xi2, 1) + (F_3y*zsc)*macaulay(i, x3, 1)) - P*(ha/2)*np.cos(theta)*macaulay(i, xi2, 1) - P*zsc*np.sin(theta)*macaulay(i, xi2, 1) + AV3.TripleIntegralZSC(i, zsc)) + c5)

    return x, y


