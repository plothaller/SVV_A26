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
		self.str_area = width_str * thickness_str + height_str * thickness_str - (math.pow(thickness_str,2))
		print("The stringer area is:", self.str_area)
		self.perimeter = math.pi*height/2 + 2*math.sqrt(math.pow(height/2,2) + math.pow(chord - height/2,2))
		self.spacing = (math.pi*height/2 + 2*math.sqrt(math.pow(height/2,2) + math.pow(chord - height/2,2)))/(number_str)
		print("The spacing is:", self.spacing)
		self.lenght_skin = np.sqrt(math.pow(self.height/2,2)+math.pow(self.chord - self.height/2,2))
		print("The skin lenght is:", self.lenght_skin)
		self.spacing_extra_booms = self.spacing/booms_per_str
		self.booms_per_str = booms_per_str
		self.a1 = math.pi*math.pow(self.height/2,2)
		self.a2 = (self.height/2) * (self.chord-self.height/2)
		self.booms_z, self.booms_y, self.booms_area = self.booms(self.booms_per_str, True)[0], self.booms(self.booms_per_str, True)[1], self.booms(self.booms_per_str, True)[2]
		self.shear_flow_magnitude = np.zeros(len(self.booms_z)+5)
		self.shear_flow_integrated = np.zeros(len(self.booms_z)+5)
		self.centroid_z = self.centroid()[0]
		self.centroid_y = self.centroid()[1]
		self.I_zz = self.moments_of_inertia(self.booms_z, self.booms_y)[0]
		print("The I_zz is:", self.I_zz)
		self.I_yy = self.moments_of_inertia(self.booms_z, self.booms_y)[1]
		print("The I_yy is:", self.I_yy)

		self.shear_center_z = self.shear_center()


	def booms(self, booms_per_str, Booms = False):
		if self.plane == "CRJ700":
			z = [-0.08650000000000001, -0.049241227826901254, 0.03039227466598985, 0.11197176918465875, 0.19355126370332768, 0.27513075822199656, 0.3567102527406655, 0.3567102527406655, 0.27513075822199656, 0.19355126370332768,0.11197176918465875,0.03039227466598985, -0.049241227826901254]  
			y = [0, 0.07111646421258025, 0.07988633519847015, 0.06213381626547678, 0.044381297332483416, 0.02662877839949005, 0.008876259466496686, -0.008876259466496686,-0.02662877839949005, -0.044381297332483416, -0.06213381626547678,-0.07988633519847015,-0.07111646421258025]
		elif self.plane == "B737":
			print("B737 still not implemented")
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
		ax.scatter(self.shear_center_z,0)
		ax.set_aspect(aspect=1)
		print("The booms are located in		X: ", self.booms_z, "	Y:",self.booms_y)
		print("The angles are for node 1:", 180/math.pi*math.acos(-self.booms_z[1]/(self.height/2)))
		plt.show()
	
	def centroid(self):
		centroid_z_area = (2*self.height)/(3*math.pi) * math.pi*self.height/2*self.skin_thickness + 2*((self.skin_thickness*self.lenght_skin) * ((self.chord-self.height/2)/2)) + 0*(self.height*self.skin_thickness)
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
		print("I_ZZ is:", I_zz)
		return I_zz, I_yy

	def base_shear_flow(self):
		self.shear_flow_magnitude = np.zeros(len(self.booms_z)+5)
		if self.plane == "CRJ700":
		#Region 1
			self.shear_flow_magnitude[0] = -1/self.I_zz*	((-math.cos(self.spacing/(self.height/2))+math.cos(0))*self.skin_thickness*math.pow(self.height/2,2))
			self.shear_flow_magnitude[1] = -1/self.I_zz*	((-math.cos(math.pi/2)+math.cos(self.spacing/(self.height/2)))* self.skin_thickness*math.pow(self.height/2,2) + 0.0711*self.str_area) + self.shear_flow_magnitude[0]
		#Region 2
			self.shear_flow_magnitude[17] = -1/self.I_zz*	(0.5*self.spar_thickness*math.pow(self.height/2,2))
		#Region 3
			c = self.spacing*(2) - math.pi*self.height/4
			self.shear_flow_magnitude[2] = -1/self.I_zz*	((self.skin_thickness*self.height *(1/2*c -1/(4*self.lenght_skin)*math.pow(c,2)))) + self.shear_flow_magnitude[1] + self.shear_flow_magnitude[17] 
			for i in range(3,7):
				c = self.spacing*(i) - math.pi*self.height/4
				self.shear_flow_magnitude[i] = -1/self.I_zz*	((self.skin_thickness * (self.height/2*c+ (-self.height/4)/(self.lenght_skin)*math.pow(c,2))) +self.sum_booms_SC(2,i-1)) + self.shear_flow_magnitude[1] + self.shear_flow_magnitude[17]
			c = self.lenght_skin
			self.shear_flow_magnitude[7] = -1/self.I_zz*	((self.skin_thickness*self.height *(1/2*c -1/(4*self.lenght_skin)*math.pow(c,2)))+self.sum_booms_SC(2,6)) + self.shear_flow_magnitude[1] + self.shear_flow_magnitude[17]
		#Region 4
			c = self.spacing*(7) - self.spacing*(6.5)
			self.shear_flow_magnitude[8] = -1/self.I_zz*	(self.skin_thickness*self.height * -1/(4*self.lenght_skin)*math.pow(c,2)) + self.shear_flow_magnitude[7]
			for i in range(9,13):
				c = self.spacing*(i-1) - self.spacing*6.5
				self.shear_flow_magnitude[i] = -1/self.I_zz*	((self.skin_thickness*self.height * -1/(4*self.lenght_skin)*math.pow(c,2))+self.sum_booms_SC(7,i-2)) + self.shear_flow_magnitude[7]
			c = self.lenght_skin
			self.shear_flow_magnitude[13] = -1/self.I_zz *	((self.skin_thickness*self.height * -1/(4*self.lenght_skin)*math.pow(c,2))+self.sum_booms_SC(7,11)) + self.shear_flow_magnitude[7]
		#Region 5
			self.shear_flow_magnitude[16] = -1/self.I_zz*	 (0.5*self.spar_thickness*math.pow(self.height/2,2))
		#Region 6
			self.shear_flow_magnitude[14] = -1/self.I_zz*	((-math.cos(-self.spacing/(self.height/2))+math.cos(-math.pi/2))*self.skin_thickness*math.pow(self.height/2,2)) + self.shear_flow_magnitude[13] - self.shear_flow_magnitude[16]
			self.shear_flow_magnitude[15] = -1/self.I_zz*	((-math.cos(0)+math.cos(-self.spacing/(self.height/2)))*self.skin_thickness*math.pow(self.height/2,2)+self.sum_booms_SC(12,12)) + self.shear_flow_magnitude[14]
			
			
		print("Base shear flow magnitude:", self.shear_flow_magnitude)
		print("Shear flow sum at the bottom joint:", self.shear_flow_magnitude[13]+self.shear_flow_magnitude[14]-self.shear_flow_magnitude[16])
		print("Shear flow sum at the top joint:", self.shear_flow_magnitude[1]+self.shear_flow_magnitude[17]-self.shear_flow_magnitude[2])
		return
	def actual_shear_flow(self):
		self.base_shear_flow()
		qs01, qs02 = self.qs0()[0], self.qs0()[1]
		print("shear flow magnitude:", self.shear_flow_magnitude)
		if self.plane == "CRJ700":
			self.shear_flow_magnitude[0] = self.shear_flow_magnitude[0] - qs01
			self.shear_flow_magnitude[1] = self.shear_flow_magnitude[1] - qs01 
			self.shear_flow_magnitude[15] = self.shear_flow_magnitude[15] - qs01
			self.shear_flow_magnitude[14] = self.shear_flow_magnitude[14] - qs01
			self.shear_flow_magnitude[16] = self.shear_flow_magnitude[16] + qs01 - qs02																				  
			self.shear_flow_magnitude[17] = self.shear_flow_magnitude[17] + qs01 - qs02
			self.shear_flow_magnitude[2] = self.shear_flow_magnitude[2] - qs02	 
			self.shear_flow_magnitude[13] = self.shear_flow_magnitude[13] - qs02
			for i in range(3,13):
				self.shear_flow_magnitude[i] = self.shear_flow_magnitude[i] - qs02
		print("After adding QS01 and QS02")
		print("Shear flow sum at the bottom joint:", self.shear_flow_magnitude[13]+self.shear_flow_magnitude[14]+self.shear_flow_magnitude[16])
		print("Shear flow sum at the top joint:", self.shear_flow_magnitude[1]+self.shear_flow_magnitude[17]-self.shear_flow_magnitude[1])
		return
	def shear_center(self):
		moment = 0
		self.actual_shear_flow()
		plate_moment_arm = (self.height/2 * (self.chord-self.height/2))/math.sqrt(math.pow(self.height/2,2)+math.pow(self.chord-self.height/2,2))
		if self.plane == "CRJ700":
			moment = self.height/2*self.spacing*(self.shear_flow_magnitude[0]+self.shear_flow_magnitude[15]) + self.height/2*(math.pi*self.height/4-self.spacing)*(self.shear_flow_magnitude[1]+self.shear_flow_magnitude[14]) + plate_moment_arm*(2*self.spacing-math.pi*self.height/4)*(self.shear_flow_magnitude[2]+self.shear_flow_magnitude[13])
			for i in range(3,13):
				if i == 7 or i == 8:
					moment = moment + self.shear_flow_magnitude[i] * plate_moment_arm * self.spacing/2
				else:
					moment = moment + self.shear_flow_magnitude[i] * plate_moment_arm * self.spacing
		shear_center_z = moment #Unit load
		print("SHEAR CENTER Z IS:", shear_center_z)
		print("shear flow magnitude:", self.shear_flow_magnitude)
		print("Printing booms locations")
		for i in range(0,len(self.booms_z)):
			print(self.booms_z[i]-self.booms_z[0], self.booms_y[i])
		return shear_center_z 

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
		#Region	1:
			self.shear_flow_integrated[0] = -1/self.I_zz*	((-math.sin(self.spacing/(self.height/2))+math.sin(0))*self.skin_thickness*math.pow(self.height/2,3)) 
			self.shear_flow_integrated[1] = -1/self.I_zz*	(self.skin_thickness*math.pow(self.height/2,3)*(-math.sin(math.pi/2)+math.sin(self.spacing/(self.height/2)))+self.sum_booms_SC(1,1)*(math.pi-self.spacing/(self.height/2)))+self.shear_flow_magnitude[0]*(math.pi-self.spacing/(self.height/2))
			self.int1 = self.shear_flow_integrated[0] + self.shear_flow_integrated[1]
		#Region 2:
			self.shear_flow_integrated[17] = 1/self.I_zz *0.5*self.spar_thickness*math.pow(self.height/2,3)/3
			self.int2 = self.shear_flow_integrated[17] 
		#Region 3:
			self.shear_flow_integrated[2] = -1/self.I_zz*	(self.skin_thickness*self.height/2*(0.5*math.pow(2*self.spacing-math.pi*self.height/4,2) - 1/6 * math.pow(2*self.spacing-math.pi*self.height/4,3)/self.lenght_skin))+(self.shear_flow_magnitude[1]+self.shear_flow_magnitude[17])*(2*self.spacing-math.pi*self.height/4)
			self.shear_flow_integrated[3] = -1/self.I_zz*	(self.skin_thickness*self.height/2*(0.5*math.pow(3*self.spacing-math.pi*self.height/4,2) - 1/6 * math.pow(3*self.spacing-math.pi*self.height/4,3)/self.lenght_skin - (0.5*math.pow(2*self.spacing-math.pi*self.height/4,2) - 1/6 * math.pow(2*self.spacing-math.pi*self.height/4,3)/self.lenght_skin)) + self.sum_booms_SC(2,2)*self.spacing) + self.shear_flow_magnitude[2]*self.spacing
			self.shear_flow_integrated[4] = -1/self.I_zz*	(self.skin_thickness*self.height/2*(0.5*math.pow(4*self.spacing-math.pi*self.height/4,2) - 1/6 * math.pow(4*self.spacing-math.pi*self.height/4,3)/self.lenght_skin - (0.5*math.pow(3*self.spacing-math.pi*self.height/4,2) - 1/6 * math.pow(3*self.spacing-math.pi*self.height/4,3)/self.lenght_skin)) + self.sum_booms_SC(3,3)*self.spacing) + self.shear_flow_magnitude[3]*self.spacing
			self.shear_flow_integrated[5] = -1/self.I_zz*	(self.skin_thickness*self.height/2*(0.5*math.pow(5*self.spacing-math.pi*self.height/4,2) - 1/6 * math.pow(5*self.spacing-math.pi*self.height/4,3)/self.lenght_skin - (0.5*math.pow(4*self.spacing-math.pi*self.height/4,2) - 1/6 * math.pow(4*self.spacing-math.pi*self.height/4,3)/self.lenght_skin)) + self.sum_booms_SC(4,4)*self.spacing) + self.shear_flow_magnitude[4]*self.spacing
			self.shear_flow_integrated[6] = -1/self.I_zz*	(self.skin_thickness*self.height/2*(0.5*math.pow(6*self.spacing-math.pi*self.height/4,2) - 1/6 * math.pow(6*self.spacing-math.pi*self.height/4,3)/self.lenght_skin - (0.5*math.pow(5*self.spacing-math.pi*self.height/4,2) - 1/6 * math.pow(5*self.spacing-math.pi*self.height/4,3)/self.lenght_skin)) + self.sum_booms_SC(5,5)*self.spacing) + self.shear_flow_magnitude[5]*self.spacing
			self.shear_flow_integrated[7] = -1/self.I_zz*	(self.skin_thickness*self.height/2*(0.5*math.pow(6.5*self.spacing-math.pi*self.height/4,2) - 1/6 * math.pow(6.5*self.spacing-math.pi*self.height/4,3)/self.lenght_skin - (0.5*math.pow(6*self.spacing-math.pi*self.height/4,2) - 1/6 * math.pow(6*self.spacing-math.pi*self.height/4,3)/self.lenght_skin)) + self.sum_booms_SC(6,6)*0.5*self.spacing) + self.shear_flow_magnitude[6]*0.5*self.spacing
			self.int3 = self.shear_flow_integrated[2] + self.shear_flow_integrated[3] + self.shear_flow_integrated[4] + self.shear_flow_integrated[5] + self.shear_flow_integrated[6] + self.shear_flow_integrated[7]

		#Region 4:
			self.shear_flow_integrated[8] = -1/self.I_zz*	(self.skin_thickness*self.height/2/self.lenght_skin*-1/6*math.pow(0.5*self.spacing,3)) + self.shear_flow_magnitude[7]*0.5*self.spacing
			self.shear_flow_integrated[9] = -1/self.I_zz*	(self.skin_thickness*self.height/2/self.lenght_skin*-1/6*(math.pow(1.5*self.spacing,3) - math.pow(0.5*self.spacing,3)) + self.sum_booms_SC(7,7))*self.spacing + self.shear_flow_magnitude[8]*self.spacing
			self.shear_flow_integrated[10] = -1/self.I_zz*	(self.skin_thickness*self.height/2/self.lenght_skin*-1/6*(math.pow(2.5*self.spacing,3) - math.pow(1.5*self.spacing,3)) + self.sum_booms_SC(8,8))*self.spacing + self.shear_flow_magnitude[9]*self.spacing
			self.shear_flow_integrated[11] = -1/self.I_zz*	(self.skin_thickness*self.height/2/self.lenght_skin*-1/6*(math.pow(3.5*self.spacing,3) - math.pow(2.5*self.spacing,3)) + self.sum_booms_SC(9,9))*self.spacing + self.shear_flow_magnitude[10]*self.spacing
			self.shear_flow_integrated[12] = -1/self.I_zz*	(self.skin_thickness*self.height/2/self.lenght_skin*-1/6*(math.pow(4.5*self.spacing,3) - math.pow(3.5*self.spacing,3)) + self.sum_booms_SC(10,10))*self.spacing + self.shear_flow_magnitude[11]*self.spacing
			self.shear_flow_integrated[13] = -1/self.I_zz*	(self.skin_thickness*self.height/2/self.lenght_skin*-1/6*(math.pow(self.lenght_skin,3) - math.pow(4.5*self.spacing,3)) + self.sum_booms_SC(11,11))*(self.lenght_skin-4.5*self.spacing) + self.shear_flow_magnitude[12]*(self.lenght_skin-4.5*self.spacing)
			self.int4 = self.shear_flow_integrated[8] + self.shear_flow_integrated[9] + self.shear_flow_integrated[10] + self.shear_flow_integrated[11] + self.shear_flow_integrated[12] + self.shear_flow_integrated[13] 
		#Region 5:
			self.shear_flow_integrated[16] = 1/self.I_zz *0.5*self.spar_thickness*math.pow(self.height/2,3)/3
			self.int5 = self.shear_flow_integrated[16]
			
		#Region 6:
			self.shear_flow_integrated[14] = -1/self.I_zz*	((-math.sin(-self.spacing/(self.height/2))+math.sin(-math.pi/2))*self.skin_thickness*math.pow(self.height/2,3)) + (self.shear_flow_magnitude[13]-self.shear_flow_magnitude[16])*(-self.spacing/(self.height/2)+math.pi/2)
			self.shear_flow_integrated[15] = -1/self.I_zz*	((-math.sin(0)+math.sin(-self.spacing/(self.height/2)))*self.skin_thickness*math.pow(self.height/2,3)+self.sum_booms_SC(12,12)*(self.spacing/(self.height/2)))+self.shear_flow_magnitude[14]*(self.spacing/(self.height/2))
			self.int6 = self.shear_flow_integrated[14] + self.shear_flow_integrated[15]

			
	def qs0(self):#in this part of the code the only thing that needs to be added is the sum of the shear flows through the arc, sum of the shear flows through the straight part of the skin and the shear flow through the spar
		self.integrate_shear_flows()
		radius_arc = self.height/2 #defining the radius of the front section
		length_straight_skin = 2*self.lenght_skin 
		sum_shearflow_through_arc = self.int1 + self.int6
		sum_shearflow_through_straightskin = self.int3 + self.int4
		shearflow_spar = self.int2 + self.int5
		
		A = np.matrix([[np.pi*self.height/2/self.skin_thickness + self.height/self.spar_thickness, - self.height/self.spar_thickness],
                               [- self.height/self.spar_thickness, length_straight_skin/self.skin_thickness + self.height/self.spar_thickness]])
		B = np.matrix([[sum_shearflow_through_arc/self.skin_thickness - shearflow_spar/self.spar_thickness],
                               [sum_shearflow_through_straightskin/self.skin_thickness + shearflow_spar/self.spar_thickness]])
		qs0 = np.linalg.solve(A,B)
		qs01 = -qs0[0]
		qs02 = -qs0[1]
		print("QS01:", qs01, "QS02:", qs02)
		return qs01, qs01


x = Geometry(17.3/100,1.1/1000,2.5/1000,1.2/1000,1.4/100,1.8/100,0.484,13,1, "CRJ700")
x.idealization()
print("Shear center location is:", x.shear_center())


#x = Geometry(20.5/100,1.1/1000,2.8/1000,1.2/1000,1.6/100,1.9/100,0.605,15,1, "B737")
#x.idealization()
#print("Shear center location is:", x.shear_center())

