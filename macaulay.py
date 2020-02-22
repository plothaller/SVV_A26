import numpy as np

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
	DOUBLEINTEGRAL = 0
	DOULBEINTEGRALPLUSZMINUSZSC = 0
	TRIPLEINTEGRALPLUSZMINUSZSC = 0
	FIVEINTEGRAL = 0

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
	A[2, 0] = -zsc
	A[2, 1] = -zsc
	A[2, 2] = -zsc
	A[2, 4] = -ha / 2 * np.cos(theta) + (ha / 2 - zsc) * np.sin(theta)
	#row 3
	A[3, 3] = np.cos(theta) * (la - xi1)
	A[3, 4] = (la - x1)
	A[3, 5] = (la - x2)
	A[3, 6] = (la - x3)
	#row 4
	A[4, 0] = (la - x1)
	A[4, 1] = (la - x2)
	A[4, 2] = (la - x3)
	A[4, 3] = np.sin(theta) * (la - xi1)
	#row 5
	A[5, 7] = x1
	A[5, 8] = 1
	A[5, 11] = zsc
	#row 6
	A[6, 9] = x1
	A[6, 10] =  1
	#row 7
	A[7, 0] = (x2 - x1) ** 3 / (6 * E * I_zz)
	A[7, 7] = x2
	A[7, 8] = 1
	A[7, 11] = -zsc
	A[7, 3] = np.sin(theta) * (x2 - xi1) ** 3 /(6 * E * I_zz) - zsc/(G * J) * (ha / 2 * np.cos(theta) * (x2 - xi1) - (ha / 2 - zsc) * np.sin(theta) * (x2 - xi1))
	#row 8
	A[8, 0] = -zsc ** 2 * (x2 - x1)
	A[8, 9] = x2
	A[8, 10] = 1
	A[8, 3] = np.cos(theta) * (x2 - xi1) ** 3 / (6 * E * I_zz)
	#row 9
	A[9, 0] = zsc ** 2 * (x3 - x1) / (G * J) - (x3 - x1) ** 3 / (6 * E * I_zz)
	A[9, 1] = zsc ** 2 * (x3 - x2) / (G * J) - (x3 - x2) ** 3 / (6 * E * I_zz)
	A[9, 3] = zsc / (G * J) * (ha / 2 * np.cos(theta) * (x3 - xi1) - (ha / 2 - zsc) * np.sin(theta) * (x3 - xi1)) - np.sin(theta) * (x3 - xi1) ** 3 / (6 * E * I_zz)
	A[9, 7] = x3
	A[9, 8] = -1
	A[9, 11] = 1
	#row 10
	A[10, 3] = (x3 - xi1) ** 3 / (6 * E * I_yy)
	A[10, 4] = (x3 - x1) ** 3 / (6 * E * I_yy)
	A[10, 5] = (x3 - x2) ** 3 / (6 * E * I_yy)
	A[10, 9] = x3
	A[10, 10] = 1
	#row 11
	A[11, 0] = -np.sin(theta) * (xi1 - x1) ** 3 / (6 * E * I_zz) + ha * np.cos(theta) * zsc * (xi1 - x1) ** 3 / (2 * G * J) + zsc ** 2 * np.sin(theta) * (xi1 - x1) / (G * J);
	A[11, 4] = np.cos(theta) * (xi1 - x1) ** 2 / (6 * E * I_yy)
	A[11, 7] = xi1
	A[11, 8] = 1
	A[11, 9] = xi1
	A[11, 10] = 1
	A[11, 11] = 2

	#b matrix
	b[0] = P * np.sin(theta) + DOUBLEINTEGRAL
	b[1] = P * np.cos(theta)
	b[2] = (-ha / 2) * P * np.cos(theta) + (ha / 2 - zsc) * P * np.sin(theta) + DOULBEINTEGRALPLUSZMINUSZSC
	b[3] = P * np.cos(theta) * (la - xi2)
	b[4] = P * np.sin(theta) * (la - xi2) + DOUBLEINTEGRAL
	b[5] = d1*np.cos(theta) + 1 / (E * I_zz) * FIVEINTEGRAL + zsc / (G * J) * TRIPLEINTEGRALPLUSZMINUSZSC
	b[6] = -d1*np.sin(theta)
	b[7] = 1 / (E * I_zz) * FIVEINTEGRAL + zsc / (G * J) * TRIPLEINTEGRALPLUSZMINUSZSC
	b[9] = d3*np.cos(theta) - zsc / (G * J) * TRIPLEINTEGRALPLUSZMINUSZSC + zsc / (G * J) * ha / 2 * P * np.cos(theta) * (x3 - xi2) - zsc / (G * J) * (ha / 2 - zsc) * P * np.sin(theta) * (x3 - xi2) - 1 / (E * I_zz) * FIVEINTEGRAL - P * np.sin(theta) * (x3 - xi2) ** 3 / (6 * E * I_zz)
	b[10] = -d3*np.sin(theta) + P * np.cos(theta) * (x3 - xi2) ** 3 / (6 * E * I_yy)
	b[11] = -np.sin(theta) / (E * I_zz) * FIVEINTEGRAL - ha * np.cos(theta) / (2 * G * J) * TRIPLEINTEGRALPLUSZMINUSZSC - zsc * np.sin(theta) / (G * J) * TRIPLEINTEGRALPLUSZMINUSZSC;

	# solve for x
	x = np.linalg.solve(A,b)
	print(x)

	Fy = x[0] + x[1] + x[2] + x[3]*np.sin(theta) - P*np.sin(theta)
	Fz = x[3]*np.cos(theta) - P*np.cos(theta) + x[4] + x[5] + x[6]

	return [Fy, Fz]