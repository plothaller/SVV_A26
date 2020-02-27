import unittest
import macaulay
import numpy as np
import deflections
import math as m
import numpy.testing as npt
import geometry_analytical as geom
import Aerodynamic_Loading_Main_V3 as AV3
#in node 1, shear flow goes to zero
#test integrals
#test moment of inertias and centroids
#test geometry stuff

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

class TestMacaulay(unittest.TestCase):


    def test_sum_forces(self):
        #Testing that equlibrium equations equate to zero
        x, A, b = macaulay.Macaulay(la, x1, x2, x3, xa, ha, d1, d3, theta, P, zsc, E, J, G, I_zz, I_yy)
        F_1y, F_2y, F_3y, F_I, F_1z, F_2z, F_3z, c1, c2, c3, c4, c5 = x[0][0][0], x[0][1][0], x[0][2][0], x[0][3][0], \
                                                                      x[0][4][0], x[0][5][0], x[0][6][0], x[0][7][0], \
                                                                      x[0][8][0], x[0][9][0], x[0][10][0], x[0][11][0]

        Fy = F_1y + F_2y + F_3y + F_I * np.sin(theta) - P * np.sin(theta)
        Fz = F_I * np.cos(theta) - P * np.cos(theta) + F_1z + F_2z + F_3z

        assert np.isclose(Fy, 0)
        assert np.isclose(Fz, 0)

    def test_singular_matrix(self):
        #Testing if the A matrix is singular
        x, A, b = macaulay.Macaulay(la, x1, x2, x3, xa, ha, d1, d3, theta, P, zsc, E, J, G, I_zz, I_yy)
        assert np.linalg.det(A) != 0

class Testdeflection(unittest.TestCase):

    def test_increase_la(self):
        #Testing that an increased aileron length, outputs higher deflections
        x, A, b= macaulay.Macaulay(la, x1, x2, x3, xa, ha, d1, d3, theta, P, zsc, E, J, G, I_zz, I_yy)
        x_1, A, b = macaulay.Macaulay(4, x1, x2, x3, xa, ha, d1, d3, theta, P, zsc, E, J, G, I_zz, I_yy)
        F_1y, F_2y, F_3y, F_I, F_1z, F_2z, F_3z, c1, c2, c3, c4, c5 = x[0][0][0], x[0][1][0], x[0][2][0], x[0][3][0], \
                                                                      x[0][4][0], x[0][5][0], x[0][6][0], x[0][7][0], \
                                                                      x[0][8][0], x[0][9][0], x[0][10][0], x[0][11][0]
        F_1y1, F_2y1, F_3y1, F_I1, F_1z1, F_2z1, F_3z1, c11, c21, c31, c41, c51 = x_1[0][0][0], x_1[0][1][0], x_1[0][2][0], x_1[0][3][0], \
                                                                      x_1[0][4][0], x_1[0][5][0], x_1[0][6][0], x_1[0][7][0], \
                                                                      x_1[0][8][0], x_1[0][9][0], x_1[0][10][0], x_1[0][11][0]

        x, y = deflections.v_deflection(la, F_I, F_1y, F_2y, F_3y, x1, x2, x3, xa, theta, P, c1, c2, E, I_zz)
        x_new, y_new = deflections.v_deflection(4, F_I1, F_1y1, F_2y1, F_3y1, x1, x2, x3, xa, theta, P, c11, c21, E, I_zz)

        assert y[-1] < y_new[-1]

    def test_hinge_boundary_conditions(self):
        #Testing to see if when no load is applied, boundary conditions are still adhered too
        x, A, b = macaulay.Macaulay(la, x1, x2, x3, xa, ha, d1, d3, theta, P, zsc, E, J, G, I_zz, I_yy)
        F_1y, F_2y, F_3y, F_I, F_1z, F_2z, F_3z, c1, c2, c3, c4, c5 = x[0][0][0], x[0][1][0], x[0][2][0], x[0][3][0], \
                                                                      x[0][4][0], x[0][5][0], x[0][6][0], x[0][7][0], \
                                                                      x[0][8][0], x[0][9][0], x[0][10][0], x[0][11][0]

        x, y = deflections.v_deflection(la, F_I, F_1y, F_2y, F_3y, x1, x2, x3, xa, theta, P, c1, c2, E, I_zz)

        i = (np.abs(x - x1)).argmin()
        k = (np.abs(x - x2)).argmin()
        j = (np.abs(x - x3)).argmin()

        assert m.isclose(y[i], d1 * np.cos(theta), rel_tol=0.01)
        assert m.isclose(y[k], 0, abs_tol=0.01)
        assert m.isclose(y[j], d3 * np.cos(theta), rel_tol=0.01)



    def test_zero_deflection(self):
        #Test to make sure a unloaded aileron produces zero deflections.
        #The boundary conditions also need to set to zero
        x, A, b = macaulay.Macaulay(la, x1, x2, x3, xa, ha, 0, 0, theta, 0, zsc, E, J, G, I_zz, I_yy)
        F_1y, F_2y, F_3y, F_I, F_1z, F_2z, F_3z, c1, c2, c3, c4, c5 = x[0][0][0], x[0][1][0], x[0][2][0], x[0][3][0], \
                                                                      x[0][4][0], x[0][5][0], x[0][6][0], x[0][7][0], \
                                                                      x[0][8][0], x[0][9][0], x[0][10][0], x[0][11][0]

        x, y = deflections.v_deflection(la, F_I, F_1y, F_2y, F_3y, x1, x2, x3, xa, theta, 0, c1, c2, E, I_zz)

        assert min(y) == max(y) == 0

class TestGeometricalProperties(unittest.TestCase):

    geometry = geom.Geometry(ha, tsk, tsp, tst, hst, wst, Ca, nst, 1)

    def test_boom_locations(self):
        #test that the locations of all booms are where they're expected to be
        (locations_z, locations_y) = ((self.obj.booms_z + ha/2), self.obj.booms_y)
        #(expected_locations_z, expected_locations_y) =

        verification_location = [[-0.,  0.],
        [-0.03725877,  0.07111646],
        [-0.11689227,  0.07988634],
        [-0.19847177,  0.06213382],
        [-0.28005126,  0.0443813 ],
        [-0.36163076,  0.02662878],
        [-0.44321025,  0.00887626],
        [-0.44321025, -0.00887626],
        [-0.36163076, -0.02662878],
        [-0.28005126, -0.0443813 ],
        [-0.19847177, -0.06213382],
        [-0.11689227, -0.07988634],
        [-0.03725877, -0.07111646]]

        print(verification_location)
# 		assert 
# 	def test_boom_areas(self):
# 		
# 	def test_centroid_location(self):
# 		
# 	def test_moment_of_intertia_booms(self):
# 		
# 	def test_moment_of_inertia_plate(self):
# 	
# 	def test_moment_of_inertia_spar(self):
# 		
# 	def test_moment_of_inertia_semicircle(self):
# 		
# 	
# class TestShear(self):
# 	
# 	def test_base_shear_centre(self):
# 		#assert that location of shear centre in [condition] is where it should be
# 		
# 		
# 	def test_shear_stresses(self):
# 		#assert that shear stresses in [condition] are what they should be 
# =============================================================================
		
class TestInterpolation_Integration(unittest.TestCase):
	def test_interpolation_1(self):
		assert abs(AV3.LinearInterpolatePos(2, 4, 1, 2, 1.5) - 3) < 0.00001
		
	def test_interpolation_2(self):
		assert abs(AV3.LinearInterpolatePos(5, 8, 4, 9, 6) - 6.2) < 0.00001
		
	def test_integration_linear_1(self):
		#integrate the function y = 1 + x sampled at 0, 1, 3 until x = 3
		assert abs(AV3.integrate_1d([0, 1, 3], [1, 2, 4], 3) - 7.5) < 0.00001
		
	def test_integration_linear_2(self):
		#integrate the function y = 1 + x sampled at 0, 1, 3 until x = 2
		assert abs(AV3.integrate_1d([0, 1, 3], [1, 2, 4], 2) - 4) < 0.00001
		
	def test_integration_cubic(self):
		#Integrate the function y = -1 + 2x + x^3 until x = 2.5
		assert abs(AV3.integrate_1d([0, 1, 2, 3, 4], [-1, 2, 11, 32, 71], 2.5) - 15.125) < 0.00001 #actual value of y integrated is 13.515625
		
	def test_double_integral(self):
		#Integrate the function y = 2 - x twice over x 
		y_new, x_new = AV3.integrate_1d_list([3, 5, 6], [-1, -3, -4], 5.2)
		assert abs(AV3.integrate_1d(x_new, y_new, 5.2) - (-4.862)) < 0.00001
		
	def test_interpolate_plot(self):
		points_x, points_y = [], []
		data_points_x, data_points_y = [0, 3, 5], [2, 4.8, 3]
		for x_i in np.arange(0, 5, 0.05):
			points_x.append(x_i)
			if x_i < 3:
				points_y.append(AV3.LinearInterpolatePos(2, 4.8, 0, 3, x_i))
			else:
				points_y.append(AV3.LinearInterpolatePos(4.8, 3, 3, 5, x_i))
		plt.figure(0)
		plt.scatter(points_x, points_y, color = 'red', s = 10)
		plt.scatter(data_points_x, data_points_y, color = 'blue', s = 60)
		plt.show()

if __name__ == '__main__':
    unittest.main()
