import matplotlib.pyplot as plt
import numpy as np
import math


'''
TODO:

-Coordinate system uses X-Y (code) instead on Z-Y (report)
-Shear center
-Torsional constant
-Centroid
-Fix boom spacing
-What happens when the boom is exactly in the vertical plate and leading edge and trailing edge plate joint?????????
-Trailing edge shear element is not considered
-Aliron height is the full spar or only half the spar?????
-Shear flow cases 2,5 (spar) are yet not implemented
-Do we want to implement skin into the booms with a large number of booms (10k-100k) ????
-Shear flows are wrong, they are just the integral component
-Shear flow calculation uses wrong integration bounds
-On the circle (L.E) to plates interface, take into accunt the lenght of the plate between the last boom in the circle and the vertical plate.
-Be careful with the inputs, the height is the full aliron thickness, not just half

'''


class Geometry:

	def __init__(self, height, skin_t, spar_t, thickness_str, height_str, width_str, chord, number_str, booms_per_str, plane):
		self.plane = plane
		self.height = height			#height of aileron
		self.skin_thickness = skin_t	#Thickness of skin
		self.spar_thickness = spar_t	#Thickness of spar
		self.t_st = thickness_str		#thickness of stringer
		self.h_st = height_str			#height of stringer
		self.w_st = width_str			#width of stringer
		self.chord = chord				#chord length
		self.n_str = number_str			#number of stringers
		self.G = 28 * math.pow(10,9)
		self.E = 73.1*math.pow(10,9)
		self.str_area = width_str * thickness_str + height_str * thickness_str - (math.pow(thickness_str,2))
		print("The stringer area is:", self.str_area)
		self.perimeter = math.pi*height/2 + 2*math.sqrt(math.pow(height/2,2) + math.pow(chord - height/2,2))
		self.spacing = (math.pi*height/2 + 2*math.sqrt(math.pow(height/2,2) + math.pow(chord - height/2,2)))/(number_str)
		print("The spacing is:", self.spacing)
		self.lenght_skin = np.sqrt(math.pow(self.height/2,2)+math.pow(self.chord - self.height/2,2))
		print("The skin lenght is:", self.lenght_skin)
		self.spacing_extra_booms = self.spacing/booms_per_str
		self.booms_per_str = booms_per_str
		self.a1 = 0.5*math.pi*math.pow(self.height/2,2)
		print("Area1 circ:", self.a1)
		self.a2 = (self.height/2) * (self.chord-self.height/2)
		print("Area2 triang:", self.a2)
		self.booms_z, self.booms_y, self.booms_area = self.booms(self.booms_per_str, True)[0], self.booms(self.booms_per_str, True)[1], self.booms(self.booms_per_str, True)[2]
		self.shear_flow_magnitude_y = np.zeros(len(self.booms_z)+5)
		self.shear_flow_integrated_y = np.zeros(len(self.booms_z)+5)
		self.shear_flow_magnitude_z = np.zeros(len(self.booms_z)+5)
		self.centroid_z = self.centroid()[0]
		self.centroid_y = self.centroid()[1]
		self.I_zz = self.moments_of_inertia(self.booms_z, self.booms_y)[0]
		print("The I_zz is:", self.I_zz)
		self.I_yy = self.moments_of_inertia(self.booms_z, self.booms_y)[1]
		print("The I_yy is:", self.I_yy)



	def booms(self, booms_per_str, Booms = False):
		if self.plane == "CRJ700":
			z = [-0.08650000000000001, -0.049241227826901254, 0.03039227466598985, 0.11197176918465875, 0.19355126370332768, 0.27513075822199656, 0.3567102527406655, 0.3567102527406655, 0.27513075822199656, 0.19355126370332768,0.11197176918465875,0.03039227466598985, -0.049241227826901254]  
			y = [0, 0.07111646421258025, 0.07988633519847015, 0.06213381626547678, 0.044381297332483416, 0.02662877839949005, 0.008876259466496686, -0.008876259466496686,-0.02662877839949005, -0.044381297332483416, -0.06213381626547678,-0.07988633519847015,-0.07111646421258025]
		elif self.plane == "B737":
			z =  [-0, -0.03692049, -0.12081074, -0.20884515, -0.29687956, -0.38491397, -0.47294838, -0.56098279, -0.56098279, -0.47294838, -0.38491397, -0.29687956, -0.20884515, -0.12081074 , -0.03692049]    
 #[-0.03692049  0.07877549]		
 #[-0.12081074  0.09876497]				
 #[-0.20884515  0.08080771]					 
 #[-0.29687956  0.06285044]					 
 #[-0.38491397  0.04489317]					 
 #[-0.47294838  0.0269359 ]					 
 #[-0.56098279  0.00897863]					
 #[-0.56098279 -0.00897863]					
 #[-0.47294838 -0.0269359 ]				
 #[-0.38491397 -0.04489317]				
 #[-0.29687956 -0.06285044]					 
 #[-0.20884515 -0.08080771]					 
 #[-0.12081074 -0.09876497]					 
 #[-0.03692049 -0.07877549]]	
		a = []
		for i in z:
			a.append(self.str_area)
		return z,y,a

	def idealization(self):
		fig = plt.figure()
		ax = fig.add_subplot(111)
		x_boom, y_boom = self.booms_z, self.booms_y
		x_circle = np.linspace(-self.height/2, 0, 50)
		y_circle = np.sqrt((-self.height/2)**2 - (x_circle)**2)
		x_plate = np.linspace(0, self.chord - self.height/2, 50)
		y_plate = self.height/2 - (self.height/2)/(self.chord-self.height/2) * x_plate
		y_vplate = np.linspace(-self.height/2, self.height/2, 2)
		x_vplate = y_vplate*0
		ax.plot(x_circle, y_circle,'b')
		ax.plot(x_circle, -y_circle,'b')
		ax.plot(x_plate, y_plate,'b')
		ax.plot(x_plate, -y_plate,'b')
		ax.plot(x_vplate, y_vplate,'b')
		ax.scatter(x_boom, y_boom)
		ax.scatter(self.centroid_z, self.centroid_y)
		self.shear_center_z = self.shear_center()
		ax.scatter(self.shear_center_z,0)
		ax.set_aspect(aspect=1)
		print("The booms are located in		X: ", self.booms_z, "	Y:",self.booms_y)
		print("The angles are for node 1:", 180/math.pi*math.acos(-self.booms_z[1]/(self.height/2)))
		plt.show()
	
	def centroid(self):
		#centroid_z_area = (2*self.height)/(3*math.pi) * math.pi*self.height/2*self.skin_thickness + 2*((self.skin_thickness*self.lenght_skin) * ((self.chord-self.height/2)/2)) + 0*(self.height*self.skin_thickness)
		centroid_z_area = - math.pi*self.height/2*self.skin_thickness*(self.height/math.pi) + 2*self.lenght_skin*self.skin_thickness*(self.chord-self.height/2)/2
		#centroid measured from middle plate: component semi-circular area: pi*r * (-2*pi/r) ; component plate: 0 ; component rear plates: 2*(self.lenght_skin*self.thickness)(self.chord-self.height/2)/2
		centroid_y_area = 0
		total_area = math.pi*self.height/2*self.skin_thickness + 2*(self.skin_thickness*self.lenght_skin) + (self.height*self.skin_thickness)
		for i in range(0,len(self.booms_z)):
			centroid_z_area = centroid_z_area + self.booms_z[i]*self.str_area
			centroid_y_area = centroid_y_area + self.booms_y[i]*self.str_area
			total_area = total_area + self.str_area
		print("Centroid Z:", centroid_z_area/total_area)
		return centroid_z_area/total_area, centroid_y_area/total_area

	def moments_of_inertia(self, z_boom, y_boom):
		I_zz = 0
		I_yy = 0
		beta = -math.acos((self.chord-self.height/2)/(self.lenght_skin))
		I_zz = 2*((math.pow(self.lenght_skin,3)*self.skin_thickness*math.pow(math.sin(beta),2))/(12) + (self.lenght_skin*self.skin_thickness)*math.pow(self.height/4,2)) + (self.skin_thickness*math.pow(self.height,3))/12 + (math.pi*math.pow(self.height/2,3)*self.skin_thickness)/2
		for z in z_boom:
			I_yy = I_yy + math.pow(abs(z-self.centroid_z),2) * self.str_area
		for y in y_boom:
			I_zz = I_zz + math.pow(abs(y-self.centroid_y),2) * self.str_area
		I_zz = 5.8159389575991465 * math.pow(10,-6)
		I_yy = 4.3632767 * math.pow(10,-5)
		print("I_ZZ is:", I_zz)
		return I_zz, I_yy

	def base_shear_flow(self):
		if self.plane == "CRJ700":
		#Vertical shear (Y-DIR)
		#------------------------
		#Region Y 1
			self.shear_flow_magnitude_y[0] = -1/self.I_zz*	((-math.cos(self.spacing/(self.height/2))+math.cos(0))*self.skin_thickness*math.pow(self.height/2,2))
			self.shear_flow_magnitude_y[1] = -1/self.I_zz*	((-math.cos(math.pi/2)+math.cos(self.spacing/(self.height/2)))* self.skin_thickness*math.pow(self.height/2,2) + 0.0711*self.str_area) + self.shear_flow_magnitude_y[0]
		#Region Y 2
			self.shear_flow_magnitude_y[17] = -1/self.I_zz*	(0.5*self.spar_thickness*math.pow(self.height/2,2))
		#Region Y 3
			c = self.spacing*(2) - math.pi*self.height/4
			self.shear_flow_magnitude_y[2] = -1/self.I_zz*	((self.skin_thickness*self.height *(1/2*c -1/(4*self.lenght_skin)*math.pow(c,2)))) + self.shear_flow_magnitude_y[1] + self.shear_flow_magnitude_y[17] 
			for i in range(3,7):
				c = self.spacing*(i) - math.pi*self.height/4
				self.shear_flow_magnitude_y[i] = -1/self.I_zz*	((self.skin_thickness * (self.height/2*c+ (-self.height/4)/(self.lenght_skin)*math.pow(c,2))) +self.sum_booms_SC(2,i-1)) + self.shear_flow_magnitude_y[1] + self.shear_flow_magnitude_y[17]
			c = self.lenght_skin
			self.shear_flow_magnitude_y[7] = -1/self.I_zz*	((self.skin_thickness*self.height *(1/2*c -1/(4*self.lenght_skin)*math.pow(c,2)))+self.sum_booms_SC(2,6)) + self.shear_flow_magnitude_y[1] + self.shear_flow_magnitude_y[17]
		#Region Y 4
			c = self.spacing*(7) - self.spacing*(6.5)
			self.shear_flow_magnitude_y[8] = -1/self.I_zz*	(self.skin_thickness*self.height * -1/(4*self.lenght_skin)*math.pow(c,2)) + self.shear_flow_magnitude_y[7]
			for i in range(9,13):
				c = self.spacing*(i-1) - self.spacing*6.5
				self.shear_flow_magnitude_y[i] = -1/self.I_zz*	((self.skin_thickness*self.height * -1/(4*self.lenght_skin)*math.pow(c,2))+self.sum_booms_SC(7,i-2)) + self.shear_flow_magnitude_y[7]
			c = self.lenght_skin
			self.shear_flow_magnitude_y[13] = -1/self.I_zz *	((self.skin_thickness*self.height * -1/(4*self.lenght_skin)*math.pow(c,2))+self.sum_booms_SC(7,11)) + self.shear_flow_magnitude_y[7]
		#Region Y 5
			self.shear_flow_magnitude_y[16] = -1/self.I_zz*	 (0.5*self.spar_thickness*math.pow(self.height/2,2))
		#Region Y 6
			self.shear_flow_magnitude_y[14] = -1/self.I_zz*	((-math.cos(-self.spacing/(self.height/2))+math.cos(-math.pi/2))*self.skin_thickness*math.pow(self.height/2,2)) + self.shear_flow_magnitude_y[13] - self.shear_flow_magnitude_y[16]
			self.shear_flow_magnitude_y[15] = -1/self.I_zz*	((-math.cos(0)+math.cos(-self.spacing/(self.height/2)))*self.skin_thickness*math.pow(self.height/2,2)+self.sum_booms_SC(12,12)) + self.shear_flow_magnitude_y[14]

		#Horizontal shear (Z-DIR)
		#------------------------
		#Region Z 1
			self.shear_flow_magnitude_z[0] = -1/self.I_zz*	((-math.cos(self.spacing/(self.height/2))+math.cos(0))*self.skin_thickness*math.pow(self.height/2,2))
			self.shear_flow_magnitude_z[1] = -1/self.I_zz*	((-math.cos(math.pi/2)+math.cos(self.spacing/(self.height/2)))* self.skin_thickness*math.pow(self.height/2,2) + 0.0711*self.str_area) + self.shear_flow_magnitude_y[0]
		#Region Z 2
			self.shear_flow_magnitude_z[17] = -1/self.I_zz*	(0.5*self.spar_thickness*math.pow(self.height/2,2))
		#Region Z 3
			c = self.spacing*(2) - math.pi*self.height/4
			self.shear_flow_magnitude_z[2] = -1/self.I_zz*	((self.skin_thickness*self.height *(1/2*c -1/(4*self.lenght_skin)*math.pow(c,2)))) + self.shear_flow_magnitude_y[1] + self.shear_flow_magnitude_y[17] 
			for i in range(3,7):
				c = self.spacing*(i) - math.pi*self.height/4
				self.shear_flow_magnitude_z[i] = -1/self.I_zz*	((self.skin_thickness * (self.height/2*c+ (-self.height/4)/(self.lenght_skin)*math.pow(c,2))) +self.sum_booms_SC(2,i-1)) + self.shear_flow_magnitude_y[1] + self.shear_flow_magnitude_y[17]
			c = self.lenght_skin
			self.shear_flow_magnitude_z[7] = -1/self.I_zz*	((self.skin_thickness*self.height *(1/2*c -1/(4*self.lenght_skin)*math.pow(c,2)))+self.sum_booms_SC(2,6)) + self.shear_flow_magnitude_y[1] + self.shear_flow_magnitude_y[17]
		#Region Z 4
			c = self.spacing*(7) - self.spacing*(6.5)
			self.shear_flow_magnitude_z[8] = -1/self.I_zz*	(self.skin_thickness*self.height * -1/(4*self.lenght_skin)*math.pow(c,2)) + self.shear_flow_magnitude_y[7]
			for i in range(9,13):
				c = self.spacing*(i-1) - self.spacing*6.5
				self.shear_flow_magnitude_z[i] = -1/self.I_zz*	((self.skin_thickness*self.height * -1/(4*self.lenght_skin)*math.pow(c,2))+self.sum_booms_SC(7,i-2)) + self.shear_flow_magnitude_y[7]
			c = self.lenght_skin
			self.shear_flow_magnitude_z[13] = -1/self.I_zz *	((self.skin_thickness*self.height * -1/(4*self.lenght_skin)*math.pow(c,2))+self.sum_booms_SC(7,11)) + self.shear_flow_magnitude_y[7]
		#Region Z 5
			self.shear_flow_magnitude_z[16] = -1/self.I_zz*	 (0.5*self.spar_thickness*math.pow(self.height/2,2))
		#Region Z 6
			self.shear_flow_magnitude_z[14] = -1/self.I_zz*	((-math.cos(-self.spacing/(self.height/2))+math.cos(-math.pi/2))*self.skin_thickness*math.pow(self.height/2,2)) + self.shear_flow_magnitude_y[13] - self.shear_flow_magnitude_y[16]
			self.shear_flow_magnitude_z[15] = -1/self.I_zz*	((-math.cos(0)+math.cos(-self.spacing/(self.height/2)))*self.skin_thickness*math.pow(self.height/2,2)+self.sum_booms_SC(12,12)) + self.shear_flow_magnitude_y[14]
		return

	def actual_shear_flow(self):
		qs01, qs02 = self.qs0()[0], self.qs0()[1]
		if self.plane == "CRJ700":
			self.shear_flow_magnitude_y[0] = self.shear_flow_magnitude_y[0] - qs01
			self.shear_flow_magnitude_y[1] = self.shear_flow_magnitude_y[1] - qs01 
			self.shear_flow_magnitude_y[15] = self.shear_flow_magnitude_y[15] - qs01
			self.shear_flow_magnitude_y[14] = self.shear_flow_magnitude_y[14] - qs01
			self.shear_flow_magnitude_y[16] = self.shear_flow_magnitude_y[16] + qs01 - qs02																				  
			self.shear_flow_magnitude_y[17] = self.shear_flow_magnitude_y[17] + qs01 - qs02
			self.shear_flow_magnitude_y[2] = self.shear_flow_magnitude_y[2] - qs02	 
			self.shear_flow_magnitude_y[13] = self.shear_flow_magnitude_y[13] - qs02
			for i in range(3,13):
				self.shear_flow_magnitude_y[i] = self.shear_flow_magnitude_y[i] - qs02
		return

	def sum_booms_SC(self, start, end):
		summation = 0
		if end > len(self.booms_y):
			raise ValueError('Sum_boom_areas. End point is greater than the number of y_booms')
		if start == end:
			print("Summing booms for:", start)
			return self.booms_y[start] * self.str_area 
		for i in range(start, end+1):
			print("Summing booms for:", i)
			summation = summation + self.booms_y[i] * self.str_area
		return summation
	
	def integrate_shear_flows(self):
		self.base_shear_flow()
		#Vertical shear (Y-DIR)
		#-----------------------
		#Region Y 1:
		self.shear_flow_integrated_y[0] = -1/self.I_zz*	((-math.sin(self.spacing/(self.height/2))+math.sin(0))*self.skin_thickness*math.pow(self.height/2,3))
		self.shear_flow_integrated_y[1] = -1/self.I_zz*	(self.skin_thickness*math.pow(self.height/2,3)*(-math.sin(math.pi/2)+math.sin(self.spacing/(self.height/2)))+self.sum_booms_SC(1,1)*(math.pi-self.spacing/(self.height/2)))+self.shear_flow_magnitude_y[0]*(math.pi-self.spacing/(self.height/2))
		self.int1 = self.shear_flow_integrated_y[0] + self.shear_flow_integrated_y[1]
		#Region Y 2:
		self.shear_flow_integrated_y[17] = 1/self.I_zz *0.5*self.spar_thickness*math.pow(self.height/2,3)/3
		self.int2 = self.shear_flow_integrated_y[17] 
		#Region Y 3:
		self.shear_flow_integrated_y[2] = -1/self.I_zz*	(self.skin_thickness*self.height/2*(0.5*math.pow(2*self.spacing-math.pi*self.height/4,2) - 1/6 * math.pow(2*self.spacing-math.pi*self.height/4,3)/self.lenght_skin))+(self.shear_flow_magnitude_y[1]+self.shear_flow_magnitude_y[17])*(2*self.spacing-math.pi*self.height/4)
		self.shear_flow_integrated_y[3] = -1/self.I_zz*	(self.skin_thickness*self.height/2*(0.5*math.pow(3*self.spacing-math.pi*self.height/4,2) - 1/6 * math.pow(3*self.spacing-math.pi*self.height/4,3)/self.lenght_skin - (0.5*math.pow(2*self.spacing-math.pi*self.height/4,2) - 1/6 * math.pow(2*self.spacing-math.pi*self.height/4,3)/self.lenght_skin)) + self.sum_booms_SC(2,2)*self.spacing) + self.shear_flow_magnitude_y[2]*self.spacing
		self.shear_flow_integrated_y[4] = -1/self.I_zz*	(self.skin_thickness*self.height/2*(0.5*math.pow(4*self.spacing-math.pi*self.height/4,2) - 1/6 * math.pow(4*self.spacing-math.pi*self.height/4,3)/self.lenght_skin - (0.5*math.pow(3*self.spacing-math.pi*self.height/4,2) - 1/6 * math.pow(3*self.spacing-math.pi*self.height/4,3)/self.lenght_skin)) + self.sum_booms_SC(3,3)*self.spacing) + self.shear_flow_magnitude_y[3]*self.spacing
		self.shear_flow_integrated_y[5] = -1/self.I_zz*	(self.skin_thickness*self.height/2*(0.5*math.pow(5*self.spacing-math.pi*self.height/4,2) - 1/6 * math.pow(5*self.spacing-math.pi*self.height/4,3)/self.lenght_skin - (0.5*math.pow(4*self.spacing-math.pi*self.height/4,2) - 1/6 * math.pow(4*self.spacing-math.pi*self.height/4,3)/self.lenght_skin)) + self.sum_booms_SC(4,4)*self.spacing) + self.shear_flow_magnitude_y[4]*self.spacing
		self.shear_flow_integrated_y[6] = -1/self.I_zz*	(self.skin_thickness*self.height/2*(0.5*math.pow(6*self.spacing-math.pi*self.height/4,2) - 1/6 * math.pow(6*self.spacing-math.pi*self.height/4,3)/self.lenght_skin - (0.5*math.pow(5*self.spacing-math.pi*self.height/4,2) - 1/6 * math.pow(5*self.spacing-math.pi*self.height/4,3)/self.lenght_skin)) + self.sum_booms_SC(5,5)*self.spacing) + self.shear_flow_magnitude_y[5]*self.spacing
		self.shear_flow_integrated_y[7] = -1/self.I_zz*	(self.skin_thickness*self.height/2*(0.5*math.pow(6.5*self.spacing-math.pi*self.height/4,2) - 1/6 * math.pow(6.5*self.spacing-math.pi*self.height/4,3)/self.lenght_skin - (0.5*math.pow(6*self.spacing-math.pi*self.height/4,2) - 1/6 * math.pow(6*self.spacing-math.pi*self.height/4,3)/self.lenght_skin)) + self.sum_booms_SC(6,6)*0.5*self.spacing) + self.shear_flow_magnitude_y[6]*0.5*self.spacing
		self.int3 = self.shear_flow_integrated_y[2] + self.shear_flow_integrated_y[3] + self.shear_flow_integrated_y[4] + self.shear_flow_integrated_y[5] + self.shear_flow_integrated_y[6] + self.shear_flow_integrated_y[7]

		#Region Y 4:
		self.shear_flow_integrated_y[8] = -1/self.I_zz*	(self.skin_thickness*self.height/2/self.lenght_skin*-1/6*math.pow(0.5*self.spacing,3)) + self.shear_flow_magnitude_y[7]*0.5*self.spacing
		self.shear_flow_integrated_y[9] = -1/self.I_zz*	(self.skin_thickness*self.height/2/self.lenght_skin*-1/6*(math.pow(1.5*self.spacing,3) - math.pow(0.5*self.spacing,3)) + self.sum_booms_SC(7,7))*self.spacing + self.shear_flow_magnitude_y[8]*self.spacing
		self.shear_flow_integrated_y[10] = -1/self.I_zz*	(self.skin_thickness*self.height/2/self.lenght_skin*-1/6*(math.pow(2.5*self.spacing,3) - math.pow(1.5*self.spacing,3)) + self.sum_booms_SC(8,8))*self.spacing + self.shear_flow_magnitude_y[9]*self.spacing
		self.shear_flow_integrated_y[11] = -1/self.I_zz*	(self.skin_thickness*self.height/2/self.lenght_skin*-1/6*(math.pow(3.5*self.spacing,3) - math.pow(2.5*self.spacing,3)) + self.sum_booms_SC(9,9))*self.spacing + self.shear_flow_magnitude_y[10]*self.spacing
		self.shear_flow_integrated_y[12] = -1/self.I_zz*	(self.skin_thickness*self.height/2/self.lenght_skin*-1/6*(math.pow(4.5*self.spacing,3) - math.pow(3.5*self.spacing,3)) + self.sum_booms_SC(10,10))*self.spacing + self.shear_flow_magnitude_y[11]*self.spacing
		self.shear_flow_integrated_y[13] = -1/self.I_zz*	(self.skin_thickness*self.height/2/self.lenght_skin*-1/6*(math.pow(self.lenght_skin,3) - math.pow(4.5*self.spacing,3)) + self.sum_booms_SC(11,11))*(self.lenght_skin-4.5*self.spacing) + self.shear_flow_magnitude_y[12]*(self.lenght_skin-4.5*self.spacing)
		self.int4 = self.shear_flow_integrated_y[8] + self.shear_flow_integrated_y[9] + self.shear_flow_integrated_y[10] + self.shear_flow_integrated_y[11] + self.shear_flow_integrated_y[12] + self.shear_flow_integrated_y[13]
		print("SF integrated 8", self.shear_flow_integrated_y[8])
		print("SF integrated 9", self.shear_flow_integrated_y[9])
		print("SF integrated 10", self.shear_flow_integrated_y[10])
		print("SF integrated 11", self.shear_flow_integrated_y[11])
		print("SF integrated 12", self.shear_flow_integrated_y[12])
		print("SF integrated 13", self.shear_flow_integrated_y[13])
		#Region Y 5:
		self.shear_flow_integrated_y[16] = 1/self.I_zz *0.5*self.spar_thickness*math.pow(self.height/2,3)/3
		self.int5 = self.shear_flow_integrated_y[16]
			
		#Region Y 6:
		self.shear_flow_integrated_y[14] = -1/self.I_zz*	((-math.sin(-self.spacing/(self.height/2))+math.sin(-math.pi/2))*self.skin_thickness*math.pow(self.height/2,3)) + (self.shear_flow_magnitude_y[13]-self.shear_flow_magnitude_y[16])*(-self.spacing/(self.height/2)+math.pi/2)
		self.shear_flow_integrated_y[15] = -1/self.I_zz*	((-math.sin(0)+math.sin(-self.spacing/(self.height/2)))*self.skin_thickness*math.pow(self.height/2,3)+self.sum_booms_SC(12,12)*(self.spacing/(self.height/2)))+self.shear_flow_magnitude_y[14]*(self.spacing/(self.height/2))
		self.int6 = self.shear_flow_integrated_y[14] + self.shear_flow_integrated_y[15]
		return


			
	def qs0(self):#in this part of the code the only thing that needs to be added is the sum of the shear flows through the arc, sum of the shear flows through the straight part of the skin and the shear flow through the spar
		self.integrate_shear_flows()
		radius_arc = self.height/2 #defining the radius of the front section
		length_straight_skin = 2*self.lenght_skin 
		sum_shearflow_through_arc = self.int1 + self.int6
		sum_shearflow_through_straightskin = self.int3 + self.int4
		shearflow_spar = self.int2 + self.int5

		deform1 = 0.5/self.a1 * 1/self.G * ((self.int1+self.int6)/self.skin_thickness+(self.int5+self.int2)/self.spar_thickness)
		deform2 = 0.5/self.a2 * 1/self.G * ((self.int3+self.int4)/self.skin_thickness-(self.int5+self.int2)/self.spar_thickness)
		
		A = np.matrix([[np.pi*self.height/2/self.skin_thickness + self.height/self.spar_thickness, - self.height/self.spar_thickness],
                               [- self.height/self.spar_thickness, length_straight_skin/self.skin_thickness + self.height/self.spar_thickness]])
		B = np.matrix([[-deform1*2*self.a1*self.G],
                               [-deform2*2*self.a2*self.G]])

		C = np.matrix([[2*self.a1, 2*self.a2, 0],
								[self.height/self.spar_thickness+self.height*math.pi*0.5/self.skin_thickness, -self.height/self.spar_thickness, -self.G*2*self.a1],
								[-self.height/self.spar_thickness, +self.height/self.spar_thickness+2*self.lenght_skin/self.skin_thickness, -self.G*2*self.a2]])

		D = np.matrix([[1],[0],[0]])

		qs0 = np.linalg.solve(A,B)
		J_matrix = np.linalg.solve(C,D)
		self.J = 1/(J_matrix[2,0]*self.G)

		self.qs01 = float(qs0[0])
		self.qs02 = float(qs0[1])
		print("Q01:", self.qs01, "QS02:", self.qs02)
		print("Torsional stiffness:", self.J)
		#self.actual_shear_flow()
		return

	def shear_center(self):
		self.qs0()
		shear_center_z = (self.int1+self.int6)*self.height/2 + (self.int3+self.int4)* (self.height/2/self.lenght_skin) * (self.chord-self.height/2) + 2*self.a1*self.qs01 + 2*self.a2*self.qs02 - self.height/2
		print(self.height/2)
		return shear_center_z 


x = Geometry(17.3/100,1.1/1000,2.5/1000,1.2/1000,1.4/100,1.8/100,0.484,13,1, "CRJ700")
x.idealization()
#x.integrate_shear_flows()
print("Shear center location is:", x.shear_center())

z =  [-0, -0.03692049, -0.12081074, -0.20884515, -0.29687956, -0.38491397, -0.47294838, -0.56098279, -0.56098279, -0.47294838, -0.38491397, -0.29687956, -0.20884515, -0.12081074 , -0.03692049]
z2 = []
z3 = [-0.08650000000000001, -0.04957951000000001, 0.03431073999999999, 0.12234515, 0.21037955999999997, 0.29841397, 0.38644838, 0.47448278999999993, 0.47448278999999993, 0.38644838, 0.29841397, 0.21037955999999997, 0.12234515, 0.03431073999999999, -0.04957951000000001]
for a in z:
	z2.append(-a-x.height/2)
print("THE LIST Z2 IS:", z2)

#x = Geometry(20.5/100,1.1/1000,2.8/1000,1.2/1000,1.6/100,1.9/100,0.605,15,1, "B737")
#x.idealization()
#print("Shear center location is:", x.shear_center())

