'''
Group A26 SVV
'''
import numpy as np
from macaulay import *
from deflections import *

aircraft = "CRJ700" # Write either A320, F100, CRJ700 or Do228 (bear in mind capitals); this is used for aerodynamic loading
Ca = 0.484 #m
la = 1.691  # m
x1 = 0.149  # m
x2 = 0.554  # m
x3 = 1.541  # m
xa = 27.2/100  # m
ha = 17.3/100  # m
tsk = 1.1/1000  # m
tsp = 2.5/1000  # m
tst = 1.2/1000  # m
hst = 1.4/100   # m
wst = 1.8/100   # m
nst = 13  # -
d1 = 0.681/100  # m
d3 = 2.030/100  # m
theta = np.radians(26)  # rad
#P = 37.9*1000  # N
P = 0
E = 73.1*10**9 #N/m2
G = 28*10**9 #N/m2

#Entering numbers from verification model
I_zz = 5.81593895759915e-06
I_yy = 4.363276766019503e-05
zsc = -0.09185594953325858
J = 0.00018782860610613963


#Calculating Reaction forces
x = Macaulay(la, x1, x2, x3, xa, ha, d1, d3, theta, P, zsc, E, J, G, I_zz, I_yy)
F_1y, F_2y, F_3y, F_I, F_1z, F_2z, F_3z, c1, c2, c3, c4, c5 = x[0][0][0], x[0][1][0], x[0][2][0], x[0][3][0], x[0][4][0], x[0][5][0], x[0][6][0], x[0][7][0], x[0][8][0], x[0][9][0], x[0][10][0], x[0][11][0]
print(F_1y)

#Plotting deflections
x_1, y_1 = v_deflection(la, F_I, F_1y, F_2y, F_3y, x1, x2, x3, xa, theta, P, c1, c2,  E, I_zz)
x_2, y_2 = w_deflection(la, F_I, F_1z, F_2z, F_3z, x1, x2, x3, xa, theta, P, c3, c4, E, I_yy)
x_3, y_3 = internal_moment_z(la, F_I, F_1y, F_2y, F_3y, x1, x2, x3, xa, theta, P)
x_4, y_4 = internal_moment_y(la, F_I, F_1z, F_2z, F_3z, x1, x2, x3, xa, theta, P)

plt.subplot(2, 2, 1)
plt.plot(x_1, y_1, 'r')
plt.title('Deflection in y')
plt.xlabel('x (m)')
plt.ylabel('v(x)')

plt.subplot(2, 2, 2)
plt.plot(x_2, y_2, 'g')
plt.title('Deflection in z')
plt.xlabel('x (m)')
plt.ylabel('w(x)')

plt.subplot(2, 2, 3)
plt.plot(x_3, y_3, 'b')
plt.title('Internal Moment around z')
plt.xlabel('x (m)')
plt.ylabel('M_z(x)')

plt.subplot(2, 2, 4)
plt.plot(x_4, y_4, 'm')
plt.title('Internal moment around y')
plt.xlabel('x (m)')
plt.ylabel('M_y(x)')
plt.show()
