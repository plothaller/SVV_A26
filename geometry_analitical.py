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
		self.str_area = width_str * thickness_str + height_str * thickness_str
		self.perimeter = math.pi*height/2 + 2*math.sqrt(math.pow(height/2,2) + math.pow(chord - height/2,2))
		self.spacing = (math.pi*height/2 + 2*math.sqrt(math.pow(height/2,2) + math.pow(chord - height/2,2)))/(number_str)
		self.lenght_skin = np.sqrt(math.pow(self.height/2,2)+math.pow(self.chord - self.height/2,2))
		self.spacing_extra_booms = self.spacing/booms_per_str
		self.booms_per_str = booms_per_str
		self.booms_z, self.booms_y, self.booms_area = self.booms(self.booms_per_str, True)[0], self.booms(self.booms_per_str, True)[1], self.booms(self.booms_per_str, True)[2]
		self.shear_flow_magnitude = np.zeros(len(self.booms_z)+6)
		self.shear_flow_magnitude[0] = 99
		self.centroid_z = self.centroid()[0]
		self.centroid_y = self.centroid()[1]
		self.I_zz = self.moments_of_inertia(self.booms_z, self.booms_y)[0]
		self.I_yy = self.moments_of_inertia(self.booms_z, self.booms_y)[1]

		self.shear_center_z = self.shear_center()
		print("The booms are located in ++++++++++:", self.booms_z, self.booms_y)


	def booms(self, booms_per_str, Booms = False):
		z = [-0.08650000000000001, -0.049241227826901254, 0.03039227466598985, 0.11197176918465875, 0.19355126370332768, 0.27513075822199656, 0.3567102527406655, 0.3567102527406655, 0.27513075822199656, 0.19355126370332768,0.11197176918465875,0.03039227466598985, -0.049241227826901254]  
		y = [0, 0.07111646421258025, 0.07988633519847015, 0.06213381626547678, 0.044381297332483416, 0.02662877839949005, 0.008876259466496686, -0.008876259466496686,-0.02662877839949005, -0.044381297332483416, -0.06213381626547678,-0.07988633519847015,-0.07111646421258025]
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
		centroid_z_area = 0
		centroid_y_area = 0
		total_area = 0
		for i in range(0,len(self.booms_z)):
			centroid_z_area = centroid_z_area + self.booms_z[i]*self.booms_area[i]
			centroid_y_area = centroid_y_area + self.booms_y[i]*self.booms_area[i]
			total_area = total_area + self.booms_area[i]
		return centroid_z_area/total_area, centroid_y_area/total_area

	def moments_of_inertia(self, z_boom, y_boom):
		I_zz = 0
		I_yy = 0
		for z in z_boom:
			I_yy = I_yy + math.pow(abs(z-self.centroid_z),2) * self.str_area
		for y in y_boom:
			I_zz = I_zz + math.pow(abs(y-self.centroid_y),2) * self.str_area 
		return I_zz, I_yy

	def base_shear_flow(self):
		self.shear_flow_magnitude = np.zeros(len(self.booms_z)+6)
		if self.plane == "CRJ700":
			self.shear_flow_magnitude[0] = -0.431/self.I_zz * self.skin_thickness*math.pow(self.height/2,2)
			self.shear_flow_magnitude[1] = -1/self.I_zz * (self.skin_thickness*math.pow(self.height/2,2) + 0.0711*self.str_area) 
			self.shear_flow_magnitude[15] = -0.431/self.I_zz * self.skin_thickness*math.pow(self.height/2,2)
			self.shear_flow_magnitude[14] = -1/self.I_zz * (self.skin_thickness*math.pow(self.height/2,2) + 0.0711*self.str_area)
			self.shear_flow_magnitude[16] = -1/self.I_zz * self.skin_thickness*math.pow(self.height/2,2)
			self.shear_flow_magnitude[17] = -1/self.I_zz * self.skin_thickness*math.pow(self.height/2,2)
			self.shear_flow_magnitude[2] = 1/self.I_zz * ((2*self.skin_thickness*math.pow(self.height/2,2)+0.0711*self.str_area)+self.skin_thickness * (-self.height/(4*self.lenght_skin)*math.pow(2*self.spacing-math.pi*self.height/4,2)))	 #Missing small skin contribution
			self.shear_flow_magnitude[13] = 1/self.I_zz * ((2*self.skin_thickness*math.pow(self.height/2,2)+0.0711*self.str_area)+self.skin_thickness * (-self.height/(4*self.lenght_skin)*math.pow(2*self.spacing-math.pi*self.height/4,2))) #Missing small skin contribution
			self.shear_flow_magnitude[7] = 1/self.I_zz *((self.skin_thickness * (-self.height/4)/(self.lenght_skin)*math.pow(self.spacing*(6.5) - math.pi*self.height/4,2))+self.sum_booms_SC(3-1,6))
			self.shear_flow_magnitude[8] = 1/self.I_zz *((self.skin_thickness * (-self.height/4)/(self.lenght_skin)*math.pow(self.spacing*(0.5),2))+self.sum_booms_SC(3-1,6))
			for i in range(3,6):
				self.shear_flow_magnitude[i] = -1/self.I_zz *((self.skin_thickness * (-self.height/4)/(self.lenght_skin)*math.pow(self.spacing*(i+1) - math.pi*self.height/4,2))+self.sum_booms_SC(3-1,i-1))
			for i in range(8,12):
				self.shear_flow_magnitude[i] = -1/self.I_zz *((self.skin_thickness * (-self.height/4)/(self.lenght_skin)*math.pow(self.spacing*(i+1) - self.spacing*(7-3) - math.pi*self.height/4,2))+self.sum_booms_SC(8-1,i-1)) + self.shear_flow_magnitude[7]	#check the sign
		if self.plane == "CRJ700":
			self.shear_flow_magnitude[0] = -0.431/self.I_zz * self.skin_thickness*math.pow(self.height/2,2)
			self.shear_flow_magnitude[1] = -1/self.I_zz * (self.skin_thickness*math.pow(self.height/2,2) + 0.0711*self.str_area) 
			self.shear_flow_magnitude[15] = -0.431/self.I_zz * self.skin_thickness*math.pow(self.height/2,2)
			self.shear_flow_magnitude[14] = -1/self.I_zz * (self.skin_thickness*math.pow(self.height/2,2) + 0.0711*self.str_area)
			self.shear_flow_magnitude[16] = -1/self.I_zz * self.skin_thickness*math.pow(self.height/2,2)
			self.shear_flow_magnitude[17] = -1/self.I_zz * self.skin_thickness*math.pow(self.height/2,2)
			self.shear_flow_magnitude[2] = 1/self.I_zz * ((2*self.skin_thickness*math.pow(self.height/2,2)+0.0711*self.str_area)+self.skin_thickness * (-self.height/(4*self.lenght_skin)*math.pow(2*self.spacing-math.pi*self.height/4,2)))	 #Missing small skin contribution
			self.shear_flow_magnitude[13] = 1/self.I_zz * ((2*self.skin_thickness*math.pow(self.height/2,2)+0.0711*self.str_area)+self.skin_thickness * (-self.height/(4*self.lenght_skin)*math.pow(2*self.spacing-math.pi*self.height/4,2))) #Missing small skin contribution
			self.shear_flow_magnitude[7] = 1/self.I_zz *((self.skin_thickness * (-self.height/4)/(self.lenght_skin)*math.pow(self.spacing*(6.5) - math.pi*self.height/4,2))+self.sum_booms_SC(3-1,6))
			self.shear_flow_magnitude[8] = 1/self.I_zz *((self.skin_thickness * (-self.height/4)/(self.lenght_skin)*math.pow(self.spacing*(0.5),2))+self.sum_booms_SC(3-1,6))
			for i in range(3,6):
				self.shear_flow_magnitude[i] = -1/self.I_zz *((self.skin_thickness * (-self.height/4)/(self.lenght_skin)*math.pow(self.spacing*(i+1) - math.pi*self.height/4,2))+self.sum_booms_SC(3-1,i-1))
			for i in range(8,12):
				self.shear_flow_magnitude[i] = -1/self.I_zz *((self.skin_thickness * (-self.height/4)/(self.lenght_skin)*math.pow(self.spacing*(i+1) - self.spacing*(7-3) - math.pi*self.height/4,2))+self.sum_booms_SC(8-1,i-1)) + self.shear_flow_magnitude[7]	#check the sign
		return
	def actual_shear_flow(self):
		self.base_shear_flow()
		qs01, qs02 = - self.qs0()[0], - self.qs0()[1]
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
			for i in range(3,12):
				self.shear_flow_magnitude[i] = self.shear_flow_magnitude[i] - qs02
		return 
	def shear_center(self):
		self.actual_shear_flow()
		plate_moment_arm = (self.height/2 * (self.chord-self.height/2))/math.sqrt(math.pow(self.height/2,2)+math.pow(self.chord-self.height/2,2))
		if self.plane == "CRJ700":
			moment = self.height/2 * (self.shear_flow_magnitude[0]+self.shear_flow_magnitude[1]+self.shear_flow_magnitude[14]+self.shear_flow_magnitude[15])
			for i in range(3,12):
				moment = moment + self.shear_flow_magnitude[i] * plate_moment_arm
		shear_center_z = moment #Unit load
		print("SHEAR CENTER Z IS:", shear_center_z)
		print("shear flow magnitude:", self.shear_flow_magnitude)
		return shear_center_z 

	def sum_booms_SC(self, start, end):
		summation = 0
		if end > len(self.booms_y):
			raise ValueError('Sum_boom_areas. End point is greater than the number of y_booms')
		for i in range(start, end):
			summation = summation + self.booms_y[i] * self.str_area
		return summation
	


	def qs0(self):#in this part of the code the only thing that needs to be added is the sum of the shear flows through the arc, sum of the shear flows through the straight part of the skin and the shear flow through the spar
		radius_arc = self.height/2 #defining the radius of the front section
		length_straight_skin = 2*self.lenght_skin
		if self.plane == "CRJ700": 
			sum_shearflow_through_arc = self.spacing*(self.shear_flow_magnitude[0]+self.shear_flow_magnitude[15]) + (math.pi*self.height/4-self.spacing)*(self.shear_flow_magnitude[1]+self.shear_flow_magnitude[14])
			sum_shearflow_through_straightskin = (2*self.spacing-math.pi*self.height/4)*(self.shear_flow_magnitude[2]+self.shear_flow_magnitude[13]) + self.spacing/2*(self.shear_flow_magnitude[7]+self.shear_flow_magnitude[8])
			for i in range(3,6):
				sum_shearflow_through_straightskin = sum_shearflow_through_straightskin	+ self.shear_flow_magnitude[i]*self.spacing
			for i in range(8,12):
				sum_shearflow_through_straightskin = sum_shearflow_through_straightskin	+ self.shear_flow_magnitude[i]*self.spacing
			shearflow_spar = (self.shear_flow_magnitude[15]+self.shear_flow_magnitude[16])*self.height/2
		elif self.plane == "B737":
			sum_shearflow_through_arc = 0
			sum_shearflow_through_straightskin = 0
			shearflow_spar = 0
		else:
			raise ValueError('Plane is neither B737 nor CRJ700')
		
		A = np.matrix([[np.pi*self.height/self.skin_thickness + self.height/self.spar_thickness, - self.height/self.spar_thickness],
                               [- self.height/self.spar_thickness, length_straight_skin/self.skin_thickness + self.height/self.spar_thickness]])
		B = np.matrix([[sum_shearflow_through_arc/self.skin_thickness - shearflow_spar/self.spar_thickness],
                               [sum_shearflow_through_straightskin/self.skin_thickness + shearflow_spar/self.spar_thickness]])
		qs0 = np.linalg.solve(A,B)
		qs01 = qs0[0]
		qs02 = qs0[1]
		return qs01, qs01



x = Geometry(17.3/100,1.1/1000,2.5/1000,1.2/1000,1.4/100,1.4/100,0.484,13,1, "CRJ700")
x.idealization()
print("Shear center location is:", x.shear_center())

