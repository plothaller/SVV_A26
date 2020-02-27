
import numpy as np
import Aerodynamic_Loading_Main_V3 as AV3

def Macaulay(la, x1, x2, x3, xa, ha, d1, d3, theta, P, zsc, E, J, G, I_zz, I_yy):

    # Solving for the reaction forces using Ax = b
    # x represents the column vector containing all unknowns
    # x = [F_1y, F_2y, F_3y, F_I, F_1z, F_2z, F_3z, c1, c2, c3, c4, c5]^T
    # The A matrix is populated with equilibrium equations and boundary conditions
    # The A matrix is populated as follows:
    # row 0: Sum of forces in y'-axis
    # row 1: Sum of forces in z'-axis
    # row 2: Torque
    # row 3: Moment around y'-axis (Internal moment equations evaluated at la)
    # row 4: Moment around z'-axis (Internal moment equations evaluated at la)
    # row 5: v-deflection in hinge 1
    # row 6: w-deflection in hinge 1
    # row 7: v-deflection in hinge 2
    # row 8: w-deflection in hinge 2
    # row 9: v-deflection in hinge 3
    # row 10: w-deflection in hinge 3
    # row 11: horizontal deflection of actuator 1 is equal to zero in global reference frame

    #Temporary variables:
    xi1 = x2 - xa / 2
    xi2 = x2 + xa / 2
    #zsc = shearcentre()

    #Initializing empty A and  vectors
    A = np.zeros((12, 12))
    b = np.zeros((12, 1))

    #intergral placeholders
    print("a")
    x, z = AV3.make_x_and_z()
    AeroLoading = AV3.MapAeroLoading(r"C:\Users\Guille\Documents\GitHub\SVV_A26\aerodynamicloadcrj700.dat")   #INPUT FILE LOCATION FOR AERO LOADING HERE
    w_bar = AV3.make_w_bar(AeroLoading)
<<<<<<< HEAD
    input("ojdlfkldkfddf")
# =============================================================================
#     x_max_double_integral_plus_minus_zsc = 0    #INPUT X_MAX FOR THE DOUBLE INTEGRAL PLUSMINUS Z_SC HERE
#     z_sc = 0        #INPUT SHEAR CENTRE LOCATION HERE
#     tau = AV3.make_tau(z_sc, AeroLoading)
#     x_max_double_integral = 0   #IPUT X_MAX FOR THE DOUBLE INTEGRAL HERE
#     x_max_three_plus_minus_zsc_1 = 0    #INPUT X_MAX FOR THE TRIPLE INTEGRAL PLUSMINUS Z_SC HERE
#     x_max_three_plus_minus_zsc_2 = 0    #INPUT X_MAX FOR THE TRIPLE INTEGRAL PLUSMINUS Z_SC HERE
#     x_max_five_1 = 0     #INPUT X_MAX FOR THE FIVEINTEGRAL HERE
#     x_max_five_2 = 0     #INPUT X_MAX FOR THE FIVEINTEGRAL HERE
#     x_max_five_3 = 0     #INPUT X_MAX FOR THE FIVEINTEGRAL HERE
#     x_max_five_4 = 0     #INPUT X_MAX FOR THE FIVEINTEGRAL HERE
#     
#     y_three_plus_minus_zsc_2, x_three_plus_minus_zsc_2 = AV3.integrate_1d_list(x, tau, x_max_three_plus_minus_zsc_1)   
# 	#intergral placeholders
# 	DOUBLEINTEGRAL = 0
# 	DOULBEINTEGRALPLUSZMINUSZSC = 0
# 	TRIPLEINTEGRALPLUSZMINUSZSC = 0
# 	FIVEINTEGRAL = 0
# 
#     y_five_2, x_five_2 = AV3.integrate_1d_list(x, w_bar, x_max_five_1)
#     y_five_3, x_five_3 = AV3.integrate_1d_list(x_five_2, y_five_2, x_max_five_2)
#     y_five_4, x_five_4 = AV3.integrate_1d_list(x_five_3, y_five_3, x_max_five_3)
#     
#     DOUBLEINTEGRAL = AV3.integrate_1d(x, w_bar, x_max_double_integral)
#     DOULBEINTEGRALPLUSZMINUSZSC = AV3.integrate_1d(x, tau, x_max_double_integral_plus_minus_zsc)
#     TRIPLEINTEGRALPLUSZMINUSZSC = AV3.integrate_1d(x_three_plus_minus_zsc_2, y_three_plus_minus_zsc_2, x_max_three_plus_minus_zsc_2)
#     FIVEINTEGRAL = AV3.integrate_1d(x_five_4, y_five_4, x_max_five_4)
# =============================================================================
=======

>>>>>>> 6caf1945c2a5261d2d5c2e317e40dbcf5989d5f2

    #A matrix:
    # row 0
    A[0, 0] = 1
    A[0, 1] = 1
    A[0, 2] = 1
    A[0, 3] = np.sin(theta)

    #row 1
    A[1, 3] = np.cos(theta)
    A[1, 4] = 1
    A[1, 5] = 1
    A[1, 6] = 1

    #row 2
    A[2, 0] = zsc
    A[2, 1] = zsc
    A[2, 2] = zsc
    A[2, 3] = ha / 2 * np.cos(theta) + zsc * np.sin(theta)



    #row 3
    A[3, 3] = np.cos(theta)*(la - xi1)
    A[3, 4] = (la - x1)
    A[3, 5] = (la - x2)
    A[3, 6] = (la - x3)


    #row 4
    A[4, 0] = (la - x1)
    A[4, 1] = (la - x2)
    A[4, 2] = (la - x3)
    A[4, 3] = np.sin(theta)*(la - xi1)


    #row 5
    A[5, 7] = x1
    A[5, 8] = 1
    A[5, 11] = zsc

    #row 6
    A[6, 9] = x1
    A[6, 10] = 1

    #row 7
    A[7, 0] = (x2 - x1)**3/(6*E*I_zz) + (zsc**2*(1/(G*J))*(x2-x1))
    A[7, 7] = x2
    A[7, 8] = 1
    A[7, 11] = +zsc
    A[7, 3] = (np.sin(theta)*(x2 - xi1)**3)/(6*E*I_zz) - (zsc/(G*J) * ((ha/2) * np.cos(theta) * (x2 - xi1))) + ((zsc)*(1/(G*J)) * np.sin(theta) * (x2 - xi1)*zsc)

    #row 8
    A[8, 4] = ((1/(6*E*I_yy))*(x2-x1)**3)
    A[8, 9] = x2
    A[8, 10] = 1
    A[8, 3] = np.cos(theta) * (x2 - xi1) ** 3 / (6 * E * I_yy)

    #row 9
    A[9, 0] = zsc ** 2 * (x3 - x1) / (G * J) + (x3 - x1) ** 3 / (6 * E * I_zz)
    A[9, 1] = zsc ** 2 * (x3 - x2) / (G * J) + (x3 - x2) ** 3 / (6 * E * I_zz)
    A[9, 3] = ((1/(6*E*I_zz))*np.sin(theta)*(x3-xi1)**3) + ((zsc)*(1/(G*J))*np.sin(theta)*(x3-xi1)*zsc) + ((ha/2)*zsc*(1/(G*J))*np.cos(theta)*(x3-xi1))
    A[9, 7] = x3
    A[9, 8] = 1
    A[9, 11] = zsc

    #row 10
    A[10, 3] = ((x3 - xi1) ** 3 / (6 * E * I_yy))*np.cos(theta)
    A[10, 4] = (x3 - x1) ** 3 / (6 * E * I_yy)
    A[10, 5] = (x3 - x2) ** 3 / (6 * E * I_yy)
    A[10, 9] = x3
    A[10, 10] = 1

    #row 11
    A[11, 0] = (zsc*(ha/2)*np.cos(theta)*(1/(G*J))*(xi1 - x1)) + ((1/(6*E*I_zz))*np.sin(theta)*(xi1-x1)**3) + (zsc**2*np.sin(theta)*(1/(G*J))*(xi1-x1))
    A[11, 4] = np.cos(theta) * (xi1 - x1) ** 3 / (6 * E * I_yy)
    A[11, 7] = (xi1*np.sin(theta))
    A[11, 8] = (np.sin(theta))
    A[11, 9] = (xi1*np.cos(theta))
    A[11, 10] = (np.cos(theta))
    A[11, 11] = ((ha/2)*np.cos(theta)) + (zsc*np.sin(theta))

    #b matrix
    b[0] = P * np.sin(theta) + AV3.DoubleIntegral(x[-1])
    b[1] = P * np.cos(theta)
    b[2] = (ha / 2) * P * np.cos(theta) + (zsc) * P * np.sin(theta) - AV3.DoubleIntegralZSC(x[-1], zsc)
    b[3] = P * np.cos(theta) * (la - xi2)
    b[4] = P * np.sin(theta) * (la - xi2) + AV3.DoubleIntegral(x[-1])
    b[5] = d1*np.cos(theta) + (np.cos(theta)/(E*I_zz))*AV3.FiveIntegral(x1) - ((zsc*np.cos(theta))/(G*J))*AV3.TripleIntegralZSC(x1, zsc)
    b[6] = -d1*np.sin(theta)
    b[7] = (1/(E*I_zz))*AV3.FiveIntegral(x2) - (zsc/(G*J))*AV3.TripleIntegralZSC(x2, zsc)
    b[9] = d3*np.cos(theta) + (P*np.sin(theta)*(1/(6*E*I_zz))*(x3-xi2)**3) - (zsc*(zsc)*P*np.sin(theta)*(x3-xi2)*(1/(G*J))) + (zsc*(ha/2)*(1/(G*J)*P*np.cos(theta)*(x3-xi2))) + (1/(E*I_zz))*AV3.FiveIntegral(x3) +(zsc*(1/(G*J))*AV3.TripleIntegralZSC(x3, zsc))
    b[10] = -d3*np.sin(theta) + P*np.cos(theta)*((x3 - xi2)**3)/(6*E*I_yy)
    b[11] = -((1/(G*J))*(ha/2)*np.cos(theta)*AV3.TripleIntegralZSC(xi1, zsc)) + ((1/(E*I_zz))*np.sin(theta)*AV3.FiveIntegral(xi1)) - (zsc*(1/(G*J))*np.sin(theta)*AV3.TripleIntegralZSC(xi1, zsc))

    # solve for x
    x = np.linalg.solve(A, b)
    print(A)
    print(b)
    print(x)
    return [x], A, b
