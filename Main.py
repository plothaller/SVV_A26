'''
Group A26 SVV
'''
import numpy as np
from macaulay import *
from deflections import *
from geometry_analytical import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


CRJ700 = Geometry(17.3/100,1.1/1000,2.5/1000,1.2/1000,1.4/100,1.8/100,0.484,13,1, "CRJ700")
B737 = Geometry(20.5/100,1.1/1000,2.8/1000,1.2/1000,1.6/100,1.9/100,0.605,15,1, "B737")

plane = CRJ700

von_misses, normal_stress, shear_stress, z_pos, y_pos = plane.Compute_section(2, 2, 2, 2, 2)

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
P = 37.9*1000  # N
#P = 0
E = 73.1*10**9 #N/m2
G = 28*10**9 #N/m2

#Retireving Geometry data
#geometry = Geometry(ha, tsk, tsp, tst, hst, wst, Ca, nst, 1)
#print("this is:",geometry.I_yy)

#Entering numbers from verification model
I_zz = 5.81593895759915e-06
I_yy = 4.363276766019503e-05
zsc = -0.087
J = 0.00018782860610613963


#Calculating Reaction forces

x, A, b = Macaulay(la, x1, x2, x3, xa, ha, d1, d3, theta, P, zsc, E, J, G, I_zz, I_yy)
F_1y, F_2y, F_3y, F_I, F_1z, F_2z, F_3z, c1, c2, c3, c4, c5 = x[0][0][0], x[0][1][0], x[0][2][0], x[0][3][0], x[0][4][0], x[0][5][0], x[0][6][0], x[0][7][0], x[0][8][0], x[0][9][0], x[0][10][0], x[0][11][0]



#Plotting deflections
x_1, y_1 = v_deflection(la, F_I, F_1y, F_2y, F_3y, x1, x2, x3, xa, theta, P, c1, c2,  E, I_zz)
x_2, y_2 = w_deflection(la, F_I, F_1z, F_2z, F_3z, x1, x2, x3, xa, theta, P, c3, c4, E, I_yy)
x_3, y_3 = internal_moment_z(la, F_I, F_1y, F_2y, F_3y, x1, x2, x3, xa, theta, P, zsc)
x_4, y_4 = internal_moment_y(la, F_I, F_1z, F_2z, F_3z, x1, x2, x3, xa, theta, P)
x_5, y_5 = shear_z(la, F_I, F_1z, F_2z, F_3z, x1, x2, x3, xa, theta, P)
x_6, y_6 = shear_y(la, F_I, F_1y, F_2y, F_3y, x1, x2, x3, xa, theta, P)
x_7, y_7 = torque(la, F_I, F_1y, F_2y, F_3y, x1, x2, x3, xa, theta, P, zsc, ha)
x_8, y_8 = Twist(la, F_I, F_1y, F_2y, F_3y, x1, x2, x3, xa, theta, P, c5, G, J, zsc, ha)


plt.subplot(3, 3, 1)
plt.plot(x_1, y_1, 'r')
plt.title('Deflection in y')
plt.xlabel('x (m)')
plt.ylabel('v(x)')

plt.subplot(3, 3, 2)
plt.plot(x_2, y_2, 'g')
plt.title('Deflection in z')
plt.xlabel('x (m)')
plt.ylabel('w(x)')

plt.subplot(3, 3, 3)
plt.plot(x_3, y_3, 'b')
plt.title('Internal Moment around z')
plt.xlabel('x (m)')
plt.ylabel('M_z(x)')

plt.subplot(3, 3, 4)
plt.plot(x_4, y_4, 'm')
plt.title('Internal moment around y')
plt.xlabel('x (m)')
plt.ylabel('M_y(x)')

plt.subplot(3, 3, 5)
plt.plot(x_5, y_5, 'y')
plt.title('Shear force in z')
plt.xlabel('x (m)')
plt.ylabel('S_z(x)')

plt.subplot(3, 3, 6)
plt.plot(x_6, y_6, 'k')
plt.title('Shear force in y')
plt.xlabel('x (m)')
plt.ylabel('S_y(x)')

plt.subplot(3, 3, 7)
plt.plot(x_7, y_7, 'c')
plt.title('Torque')
plt.xlabel('x (m)')
plt.ylabel('T(x)')

plt.subplot(3, 3, 8)
plt.plot(x_8, y_8, 'c')
plt.title('Twist')
plt.xlabel('x (m)')
plt.ylabel('T(x)')
plt.show()

print("DOne with macaulay")



x_pos_section = 0.5
print("The index we are using is:", x_pos_section*500/la)
print("aa", y_3[round(x_pos_section*500/la)], y_4[round(x_pos_section*500/la)], y_5[round(x_pos_section*500/la)], y_6[round(x_pos_section*500/la)], y_7[round(x_pos_section*500/la)])
#x_pos_section = 1
#print("The index we are using is:", x_pos_section*500/la)
#print("aa", y_3[round(x_pos_section*500/la)], y_4[round(x_pos_section*500/la)], y_5[round(x_pos_section*500/la)], y_6[round(x_pos_section*500/la)], y_7[round(x_pos_section*500/la)])
#Compute at section x=0.5m (to compare with the verification model)
von_misses, normal_stress, shear_stress, z_pos, y_pos = plane.Compute_section(y_3[round(x_pos_section*500/la)], y_4[round(x_pos_section*500/la)], y_5[round(x_pos_section*500/la)], y_6[round(x_pos_section*500/la)], y_7[round(x_pos_section*500/la)])
von_misses_max = np.max(von_misses)

fig = plt.figure()
ax = fig.add_subplot(111)
img = ax.scatter(z_pos, y_pos, c=von_misses, cmap=plt.cm.jet) 
fig.colorbar(img)
ax.set_aspect(aspect=1)
plt.show()



#fig = plt.figure()
#color_von_misses = plt.cm.jet(abs(von_misses/von_misses_max))
#ax = fig.add_subplot(111)
#ax.plot(z_pos, y_pos, color=color_von_misses) 
##fig.colorbar(img)
#ax.set_aspect(aspect=1)
#plt.show()

slices = 500
sections = np.linspace(0,la,slices)
return_lenght = len(plane.Compute_section(0, 0, 0, 0, 0)[0])
von_misses_total	= np.zeros(slices*return_lenght)
normal_stress_total	= np.zeros(slices*return_lenght)
shear_stress_total	= np.zeros(slices*return_lenght)
z_pos_total	= np.zeros(slices*return_lenght)
y_pos_total	= np.zeros(slices*return_lenght)
x_pos_total	= np.zeros(slices*return_lenght)
von_misses_max_total = np.zeros(slices)



for section in range(0,slices):
	print("Doing section:	[", section,"/", len(x),"]")
	von_misses, normal_stress, shear_stress, z_pos, y_pos = plane.Compute_section(y_3[section], y_4[section], y_5[section], y_6[section], y_7[section])
	von_misses_max_total[section] = np.max(von_misses) 
	for i in range(0,len(von_misses)):
		von_misses_total[section*return_lenght+i] = von_misses[i]
		z_pos_total[section*return_lenght+i] = z_pos[i]
		y_pos_total[section*return_lenght+i] = y_pos[i]
		x_pos_total[section*return_lenght+i] = section*la/slices

von_misses_max = np.max(von_misses_total)
von_misses_min = np.min(von_misses_total)
print("The maximum von misses stress is:", von_misses_max)



fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(sections[::20], von_misses_max_total[::20])

plt.show()

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
##x_boom, y_boom = self.booms_z, self.booms_y
##ax.scatter(x_boom, y_boom)
##ax.scatter(self.centroid_z, self.centroid_y)
##self.shear_center_z = self.shear_center()
##ax.scatter(self.shear_center_z,0)
img = ax.scatter(x_pos_total[::5], y_pos_total[::5]+0.15, z_pos_total[::5] ,c=von_misses_total[::5], cmap=plt.cm.RdYlGn) #cmap = ((255*(von_misses[i]-von_mises_min)/(von_mises_max-von_mises_min)),(255*(1-(von_misses[i]-von_mises_min)/(von_mises_max-von_mises_min))),0))
fig.colorbar(img)
ax.set_xlim3d(0,2)
ax.set_ylim3d(0,0.2)
ax.set_zlim3d(0,0.5)

plt.show()





