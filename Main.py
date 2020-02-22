'''
Group A26 SVV
'''
import numpy as np
from Step_function import *

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
wst = 1.8/1000   # m
nst = 13  # -
d1 = 0.681/100  # m
d3 = 2.030/100  # m
theta = np.radians(28)  # rad
P = 37.9*1000  # N

#Entering numbers from verification model
I_zz = 5.1138606931009e-06
I_yy = 3.78947094284179e-05
zsc = -0.09267209562108152
J = 0.00018782860610613963
E = 72.9*10**9
G = 27.1*10**9

Fy, Fz = Macaulay(Ca, la, x1, x2, x3, xa, ha, tsk, tsp, tst, hst, wst, nst, d1, d3, theta, P, zsc, E, J, G, I_zz, I_yy)

print(Fy)
print(Fz)