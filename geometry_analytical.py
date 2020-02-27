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
-On the circle (L.E) to plates interface, take into accunt the length of the plate between the last boom in the circle and the vertical plate.
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
		self.perimeter = math.pi*height/2 + 2*math.sqrt(math.pow(height/2,2) + math.pow(chord - height/2,2))
		self.spacing = (math.pi*height/2 + 2*math.sqrt(math.pow(height/2,2) + math.pow(chord - height/2,2)))/(number_str)
		self.length_skin = np.sqrt(math.pow(self.height/2,2)+math.pow(self.chord - self.height/2,2))
		self.spacing_extra_booms = self.spacing/booms_per_str
		self.booms_per_str = booms_per_str
		self.a1 = 0.5*math.pi*math.pow(self.height/2,2)
		self.a2 = (self.height/2) * (self.chord-self.height/2)
		self.booms_z, self.booms_y, self.booms_area = self.booms(self.booms_per_str, True)[0], self.booms(self.booms_per_str, True)[1], self.booms(self.booms_per_str, True)[2]
		self.shear_flow_magnitude_y = np.zeros(len(self.booms_z)+5)
		self.shear_flow_integrated_y = np.zeros(len(self.booms_z)+5)
		self.shear_flow_magnitude_z = np.zeros(len(self.booms_z)+5)
		self.shear_stress = np.zeros(len(self.booms_z)+5)
		self.centroid_z = self.centroid()[0]
		self.centroid_y = self.centroid()[1]
		self.I_zz = self.moments_of_inertia()[0]
		self.I_yy = self.moments_of_inertia()[1]
		if self.plane == "CRJ700":
			self.plane_delta = 0
		elif self.plane == "B737":
			self.plane_delta = 1
		else:
			raise ValueError("Plane is neither CRJ700 nor B737")



	def booms(self, booms_per_str, Booms = False):
		if self.plane == "CRJ700":
			z = [-0.08650000000000001, -0.049241227826901254, 0.03039227466598985, 0.11197176918465875, 0.19355126370332768, 0.27513075822199656, 0.3567102527406655, 0.3567102527406655, 0.27513075822199656, 0.19355126370332768,0.11197176918465875,0.03039227466598985, -0.049241227826901254]  
			y = [0, 0.07111646421258025, 0.07988633519847015, 0.06213381626547678, 0.044381297332483416, 0.02662877839949005, 0.008876259466496686, -0.008876259466496686,-0.02662877839949005, -0.044381297332483416, -0.06213381626547678,-0.07988633519847015,-0.07111646421258025]
		elif self.plane == "B737":
			z = [-0.1025, -0.06557951, 0.018310740000000006, 0.10634515000000001, 0.19437956, 0.28241397, 0.37044838, 0.45848279, 0.45848279, 0.37044838, 0.28241397, 0.19437956, 0.10634515000000001, 0.018310740000000006, -0.06557951]   
			y = [0, 0.07877549,  0.09876497, 0.08080771, 0.06285044, 0.04489317, 0.0269359, 0.00897863, -0.00897863, -0.0269359, -0.04489317, -0.06285044, -0.08080771, -0.09876497, -0.07877549]
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

	#def idealization(self):
	#	fig = plt.figure()
	#	ax = fig.add_subplot(111)
	#	x_boom, y_boom = self.booms_z, self.booms_y
	#	x_circle = np.linspace(-self.height/2, 0, 50)
	#	y_circle = np.sqrt((-self.height/2)**2 - (x_circle)**2)
	#	x_plate = np.linspace(0, self.chord - self.height/2, 50)
	#	y_plate = self.height/2 - (self.height/2)/(self.chord-self.height/2) * x_plate
	#	y_vplate = np.linspace(-self.height/2, self.height/2, 2)
	#	x_vplate = y_vplate*0
	#	ax.plot(x_circle, y_circle,'b')
	#	ax.plot(x_circle, -y_circle,'b')
	#	ax.plot(x_plate, y_plate,'b')
	#	ax.plot(x_plate, -y_plate,'b')
	#	ax.plot(x_vplate, y_vplate,'b')
	#	ax.scatter(x_boom, y_boom)
	#	ax.scatter(self.centroid_z, self.centroid_y)
	#	#self.shear_center_z = self.shear_center()
	#	#ax.scatter(self.shear_center_z,0)
	#	ax.set_aspect(aspect=1)
	#	plt.show()
	
	def centroid(self):
		#centroid_z_area = (2*self.height)/(3*math.pi) * math.pi*self.height/2*self.skin_thickness + 2*((self.skin_thickness*self.length_skin) * ((self.chord-self.height/2)/2)) + 0*(self.height*self.skin_thickness)
		#centroid_z_area = math.pi*self.height/2*self.skin_thickness*(-self.height/math.pi) + 2*self.length_skin*self.skin_thickness*(self.chord-self.height/2)/2
		total_area = math.pi*self.height/2*self.skin_thickness + self.height*self.spar_thickness + 2*self.length_skin*self.skin_thickness
		centroid_z_semicircle = math.pi*self.height/2*self.skin_thickness * (-self.height/math.pi)
		centroid_z_spar = self.height*self.spar_thickness * 0
		centroid_z_plate = self.length_skin*self.skin_thickness * (self.chord-self.height/2)/2
		centroid_z_area = centroid_z_semicircle + centroid_z_spar + 2*centroid_z_plate
	
		
		#centroid measured from middle plate: component semi-circular area: pi*r * (-2*pi/r) ; component plate: 0 ; component rear plates: 2*(self.length_skin*self.thickness)(self.chord-self.height/2)/2
		centroid_y_area = 0
		for i in range(0,len(self.booms_z)):
			centroid_z_area += self.booms_z[i]*self.str_area
			centroid_y_area += self.booms_y[i]*self.str_area
			total_area = total_area + self.str_area
		return centroid_z_area/total_area, centroid_y_area/total_area

	def moments_of_inertia(self):
		z_boom = self.booms_z
		y_boom = self.booms_y
		beta = math.acos((self.chord-self.height/2)/(self.length_skin))
		I_zz = 2*((math.pow(self.length_skin,3)*self.skin_thickness*math.pow(math.sin(beta),2))/12 + (self.length_skin*self.skin_thickness)*math.pow(self.height/4,2)) + (self.spar_thickness*math.pow(self.height,3))/12 + (math.pi/8 * ((self.height/2 + self.skin_thickness/2)**4 - (self.height/2 - self.skin_thickness/2)**4)) 	#(math.pi*math.pow(self.height/2,3)*self.skin_thickness)/2 #; alternative method for calculating the MoI of the semi-circular arc
		#I_yy = (math.pi/8*math.pow(self.height,3)*self.skin_thickness)/2 + (math.pi*self.height/2*self.skin_thickness*(self.centroid_z+(self.height/math.pi))**2) + ((self.height*self.spar_thickness**3)/12) + (self.height*self.spar_thickness*self.centroid_z**2) + 2*((math.pow(self.length_skin,3)*self.skin_thickness*math.pow(math.sin(beta),2))/(12) + (self.length_skin*self.skin_thickness) * ((self.chord-self.height/2)/2 -self.centroid_z)**2)
		#I_yy =  + 				# 2*((math.pow(self.length_skin,3)*self.skin_thickness*math.pow(math.sin(beta),2))/(12) + (self.length_skin*self.skin_thickness) * ((self.chord-self.height/2)/2 -self.centroid_z)**2)
		I_yy_semi_circular = ((math.pi/8 - 8/(9*math.pi))*((self.height/2+self.skin_thickness/2)**4 - (self.height/2 - self.skin_thickness/2)**4)) + (self.height/2*math.pi*self.skin_thickness*(self.centroid_z+self.height/math.pi)**2)
		I_yy_spar = ((self.height*self.spar_thickness**3)/12) + (self.height*self.spar_thickness*self.centroid_z**2)
		I_yy_plates = 2*(self.skin_thickness*(self.length_skin)**3 * math.cos(beta)**2/12) + ((self.chord-self.height/2)/2 - self.centroid_z)**2 * self.skin_thickness*self.length_skin*2
		I_yy = I_yy_semi_circular + I_yy_spar + I_yy_plates
		for z in z_boom:
			I_yy = I_yy + math.pow(abs(z-self.centroid_z),2) * self.str_area
		for y in y_boom:
			I_zz = I_zz + math.pow(abs(y-self.centroid_y),2) * self.str_area
		#I_zz = 5.8159389575991465 * math.pow(10,-6)
		print(" I_zz is: \t\t", I_zz, "\n reference I_zz: \t 5.8159389575991465e-06 \n percentage difference; ", (I_zz - 5.8159389575991465e-06)/(5.8159389575991465e-08), "%")
		print(" I_yy is: \t\t", I_yy, "\n reference I_yy: \t 4.363276766019503e-05 \n percentage difference; ", (I_yy - 4.363276766019503e-05)/(4.363276766019503e-07), "%")
		return I_zz, I_yy

	def base_shear_flow(self):

	#Vertical shear (Y-DIR)
	#------------------------
	#Region Y 1
		self.shear_flow_magnitude_y[0] = -1/self.I_zz*	((-math.cos(self.spacing/(self.height/2))+math.cos(0))*self.skin_thickness*math.pow(self.height/2,2))
		self.shear_flow_magnitude_y[1] = -1/self.I_zz*	((-math.cos(math.pi/2)+math.cos(self.spacing/(self.height/2)))* self.skin_thickness*math.pow(self.height/2,2) + self.booms_y[1]*self.str_area) + self.shear_flow_magnitude_y[0]
	#Region Y 2
		self.shear_flow_magnitude_y[17+2*self.plane_delta] = -1/self.I_zz*	(0.5*self.spar_thickness*math.pow(self.height/2,2))
	#Region Y 3
		c = self.spacing*(2) - math.pi*self.height/4
		self.shear_flow_magnitude_y[2] = -1/self.I_zz*	((self.skin_thickness*self.height *(1/2*c -1/(4*self.length_skin)*math.pow(c,2)))) + self.shear_flow_magnitude_y[1] + self.shear_flow_magnitude_y[17+2*self.plane_delta] 
		for i in range(3,7+self.plane_delta):
			c = self.spacing*(i) - math.pi*self.height/4
			self.shear_flow_magnitude_y[i] = -1/self.I_zz*	((self.skin_thickness * (self.height/2*c+ (-self.height/4)/(self.length_skin)*math.pow(c,2))) +self.sum_boom_SC_y(2,i-1)) + self.shear_flow_magnitude_y[1] + self.shear_flow_magnitude_y[17+2*self.plane_delta]
		c = self.length_skin
		self.shear_flow_magnitude_y[7+self.plane_delta] = -1/self.I_zz*	((self.skin_thickness*self.height *(1/2*c -1/(4*self.length_skin)*math.pow(c,2)))+self.sum_boom_SC_y(2,6+self.plane_delta)) + self.shear_flow_magnitude_y[1] + self.shear_flow_magnitude_y[17+2*+self.plane_delta]
	#Region Y 4
		c = self.spacing*0.5
		self.shear_flow_magnitude_y[8+self.plane_delta] = -1/self.I_zz*	(self.skin_thickness*self.height * -1/(4*self.length_skin)*math.pow(c,2)) + self.shear_flow_magnitude_y[7+self.plane_delta]
		for i in range(9+self.plane_delta,13+2*self.plane_delta):
			c = self.spacing*(i-1) - self.spacing*(6.5+self.plane_delta)
			self.shear_flow_magnitude_y[i] = -1/self.I_zz*	((self.skin_thickness*self.height * -1/(4*self.length_skin)*math.pow(c,2))+self.sum_boom_SC_y(7+self.plane_delta,i-2)) + self.shear_flow_magnitude_y[7+self.plane_delta]
		c = self.length_skin
		self.shear_flow_magnitude_y[13+2*self.plane_delta] = -1/self.I_zz*	((self.skin_thickness*self.height * -1/(4*self.length_skin)*math.pow(c,2))+self.sum_boom_SC_y(7+self.plane_delta,11+2*self.plane_delta)) + self.shear_flow_magnitude_y[7+self.plane_delta]
	#Region Y 5
		self.shear_flow_magnitude_y[16+2*self.plane_delta] = -1/self.I_zz*	 (0.5*self.spar_thickness*math.pow(self.height/2,2))
	#Region Y 6
		self.shear_flow_magnitude_y[14+2*self.plane_delta] = -1/self.I_zz*	((-math.cos(-self.spacing/(self.height/2))+math.cos(-math.pi/2))*self.skin_thickness*math.pow(self.height/2,2)) + self.shear_flow_magnitude_y[13+2*self.plane_delta] - self.shear_flow_magnitude_y[16+2*self.plane_delta]
		self.shear_flow_magnitude_y[15+2*self.plane_delta] = -1/self.I_zz*	((-math.cos(0)+math.cos(-self.spacing/(self.height/2)))*self.skin_thickness*math.pow(self.height/2,2)+self.sum_boom_SC_y(12+2*self.plane_delta,12+2*self.plane_delta)) + self.shear_flow_magnitude_y[14+2*self.plane_delta]
	
	#Horizontal shear (Z-DIR)
	#------------------------
	#Region Z 1
		self.shear_flow_magnitude_z[0] = -1/self.I_yy*	(self.skin_thickness*((self.height/2*math.sin(self.spacing/(self.height/2))-self.centroid_z*(self.spacing/(self.height/2)))-(self.height/2*math.sin(0)-self.centroid_z*(0)))+self.sum_boom_SC_z(0,0)*self.str_area/2)
		self.shear_flow_magnitude_z[1] = -1/self.I_yy*	(self.skin_thickness*((self.height/2*math.sin(math.pi/2)-self.centroid_z*(math.pi/2))-(self.height/2*math.sin(self.spacing/(self.height/2))-self.centroid_z*(self.spacing/(self.height/2))))+(self.booms_z[1]+self.centroid_z)*self.str_area) + self.shear_flow_magnitude_z[0] 	
	#Region Z 2
		self.shear_flow_magnitude_z[17+2*self.plane_delta] = -1/self.I_yy*	(-self.spar_thickness*self.centroid_z*self.height/2)
	#Region Z 3
		self.shear_flow_magnitude_z[2] = -1/self.I_yy*	(self.skin_thickness*((self.chord-self.height/2)/(self.length_skin)*math.pow(2*self.spacing-math.pi*self.height/4,2)-self.centroid_z*(2*self.spacing-math.pi*self.height/4))) + self.shear_flow_magnitude_z[1] + self.shear_flow_magnitude_z[17+2*self.plane_delta]
		for i in range(3,7+self.plane_delta):
			c = self.spacing*(i) - math.pi*self.height/4
			self.shear_flow_magnitude_z[i] = -1/self.I_yy* (self.skin_thickness*(((self.chord-self.height/2)/(2*self.length_skin))*math.pow(c,2)-(self.centroid_z)*c)+self.sum_boom_SC_z(2,i-1)*self.str_area) + self.shear_flow_magnitude_z[1] + self.shear_flow_magnitude_z[17+2*self.plane_delta]
		c = self.length_skin
		self.shear_flow_magnitude_z[7+self.plane_delta] = -1/self.I_yy* (self.skin_thickness*(((self.chord-self.height/2)/(2*self.length_skin))*math.pow(c,2)-(self.centroid_z)*c)+self.sum_boom_SC_z(2,6+self.plane_delta)*self.str_area) + self.shear_flow_magnitude_z[1] + self.shear_flow_magnitude_z[17+2*self.plane_delta]
	#Region Z 4
		c = self.spacing*0.5
		self.shear_flow_magnitude_z[8+self.plane_delta] = -1/self.I_yy*	(self.skin_thickness*(-(self.chord-self.height/2)/(2*self.length_skin)*math.pow(c,2)+c*(self.chord-self.height/2-self.centroid_z))) + self.shear_flow_magnitude_z[7+self.plane_delta]
		for i in range(9+self.plane_delta,13+2*self.plane_delta):
			c = self.spacing*(i-1) - self.spacing*(6.5+self.plane_delta)
			self.shear_flow_magnitude_z[i] = -1/self.I_yy*	(self.skin_thickness*(-(self.chord-self.height/2)/(2*self.length_skin)*math.pow(c,2)+c*(self.chord-self.height/2-self.centroid_z))+self.sum_boom_SC_z(7+self.plane_delta,i-2)*self.str_area) + self.shear_flow_magnitude_z[7+self.plane_delta]
		c = self.length_skin
		self.shear_flow_magnitude_z[13+2*self.plane_delta] = -1/self.I_yy *	(self.skin_thickness*(-(self.chord-self.height/2)/(2*self.length_skin)*math.pow(c,2)+c*(self.chord-self.height/2-self.centroid_z))+self.sum_boom_SC_z(7+self.plane_delta,11+2*self.plane_delta)*self.str_area) + self.shear_flow_magnitude_z[7+self.plane_delta]
	#Region Z 5
		self.shear_flow_magnitude_z[16+2*self.plane_delta] = -1/self.I_yy*	(-self.spar_thickness*self.centroid_z*self.height/2)
	#Region Z 6
		self.shear_flow_magnitude_z[14+2*self.plane_delta] = -1/self.I_yy * (self.skin_thickness*self.height/2*((-self.height/2*math.sin(-math.pi/2)-self.centroid_z*(-math.pi/2))-(-self.height/2*math.sin(-self.spacing/(self.height/2))-self.centroid_z*(-self.spacing/(self.height/2))))) + self.shear_flow_magnitude_z[13+2*self.plane_delta] - self.shear_flow_magnitude_z[16+2*self.plane_delta]
		self.shear_flow_magnitude_z[15+2*self.plane_delta] = -1/self.I_yy * (self.skin_thickness*self.height/2*((-self.height/2*math.sin(-self.spacing/(self.height/2))-self.centroid_z*(-self.spacing/(self.height/2)))-(-self.height/2*math.sin(0)-self.centroid_z*(0)))+self.sum_boom_SC_z(12+2*self.plane_delta,12+2*self.plane_delta)*self.str_area + self.sum_boom_SC_z(0,0)*self.str_area/2) + self.shear_flow_magnitude_z[14+2*self.plane_delta]
		return

	def actual_shear_flow(self):
		self.qs0()
		self.shear_center()
		if self.plane == "CRJ700":
			self.shear_flow_magnitude_y[0] = self.shear_flow_magnitude_y[0] - self.qs01
			self.shear_flow_magnitude_y[1] = self.shear_flow_magnitude_y[1] - self.qs01 
			self.shear_flow_magnitude_y[15] = self.shear_flow_magnitude_y[15] - self.qs01
			self.shear_flow_magnitude_y[14] = self.shear_flow_magnitude_y[14] - self.qs01
			self.shear_flow_magnitude_y[16] = self.shear_flow_magnitude_y[16] + self.qs01 - self.qs02																				  
			self.shear_flow_magnitude_y[17] = self.shear_flow_magnitude_y[17] + self.qs01 - self.qs02
			self.shear_flow_magnitude_y[2] = self.shear_flow_magnitude_y[2] - self.qs02	 
			self.shear_flow_magnitude_y[13] = self.shear_flow_magnitude_y[13] - self.qs02
			for i in range(3,13):
				self.shear_flow_magnitude_y[i] = self.shear_flow_magnitude_y[i] - self.qs02
		return

	def sum_boom_SC_y(self, start, end):
		summation = 0
		if end > len(self.booms_y):
			raise ValueError('Sum_boom_areas. End point is greater than the number of y_booms')
		if start == end:
			return self.booms_y[start] * self.str_area 
		for i in range(start, end+1):
			summation = summation + self.booms_y[i] * self.str_area
		return summation

	def sum_boom_SC_z(self, start, end):
		summation = 0
		if end > len(self.booms_z):
			raise ValueError('Sum_boom_areas. End point is greater than the number of y_booms')
		if start == end:
			return (self.booms_z[start]-self.centroid_z) * self.str_area 
		for i in range(start, end+1):
			summation = summation + (self.booms_z[i]-self.centroid_z) * self.str_area
		return summation
	
	def integrate_shear_flows(self):
		self.base_shear_flow()
		#Vertical shear (Y-DIR)
		#-----------------------
		#Region Y 1:
		self.shear_flow_integrated_y[0] = -1/self.I_zz*	((-math.sin(self.spacing/(self.height/2))+math.sin(0))*self.skin_thickness*math.pow(self.height/2,3))
		self.shear_flow_integrated_y[1] = -1/self.I_zz*	(self.skin_thickness*math.pow(self.height/2,3)*(-math.sin(math.pi/2)+math.sin(self.spacing/(self.height/2)))+self.sum_boom_SC_y(1,1)*(math.pi-self.spacing/(self.height/2)))+self.shear_flow_magnitude_y[0]*(math.pi-self.spacing/(self.height/2))
		self.int1 = self.shear_flow_integrated_y[0] + self.shear_flow_integrated_y[1]
		#Region Y 2:
		self.shear_flow_integrated_y[17+2*self.plane_delta] = 1/self.I_zz *0.5*self.spar_thickness*math.pow(self.height/2,3)/3
		self.int2 = self.shear_flow_integrated_y[17] 
		#Region Y 3:
		self.shear_flow_integrated_y[2] = -1/self.I_zz*	(self.skin_thickness*self.height/2*(0.5*math.pow(2*self.spacing-math.pi*self.height/4,2) - 1/6 * math.pow(2*self.spacing-math.pi*self.height/4,3)/self.length_skin))+(self.shear_flow_magnitude_y[1]+self.shear_flow_magnitude_y[17])*(2*self.spacing-math.pi*self.height/4)
		self.shear_flow_integrated_y[3] = -1/self.I_zz*	(self.skin_thickness*self.height/2*(0.5*math.pow(3*self.spacing-math.pi*self.height/4,2) - 1/6 * math.pow(3*self.spacing-math.pi*self.height/4,3)/self.length_skin - (0.5*math.pow(2*self.spacing-math.pi*self.height/4,2) - 1/6 * math.pow(2*self.spacing-math.pi*self.height/4,3)/self.length_skin)) + self.sum_boom_SC_y(2,2)*self.spacing) + self.shear_flow_magnitude_y[2]*self.spacing
		self.shear_flow_integrated_y[4] = -1/self.I_zz*	(self.skin_thickness*self.height/2*(0.5*math.pow(4*self.spacing-math.pi*self.height/4,2) - 1/6 * math.pow(4*self.spacing-math.pi*self.height/4,3)/self.length_skin - (0.5*math.pow(3*self.spacing-math.pi*self.height/4,2) - 1/6 * math.pow(3*self.spacing-math.pi*self.height/4,3)/self.length_skin)) + self.sum_boom_SC_y(3,3)*self.spacing) + self.shear_flow_magnitude_y[3]*self.spacing
		self.shear_flow_integrated_y[5] = -1/self.I_zz*	(self.skin_thickness*self.height/2*(0.5*math.pow(5*self.spacing-math.pi*self.height/4,2) - 1/6 * math.pow(5*self.spacing-math.pi*self.height/4,3)/self.length_skin - (0.5*math.pow(4*self.spacing-math.pi*self.height/4,2) - 1/6 * math.pow(4*self.spacing-math.pi*self.height/4,3)/self.length_skin)) + self.sum_boom_SC_y(4,4)*self.spacing) + self.shear_flow_magnitude_y[4]*self.spacing
		self.shear_flow_integrated_y[6] = -1/self.I_zz*	(self.skin_thickness*self.height/2*(0.5*math.pow(6*self.spacing-math.pi*self.height/4,2) - 1/6 * math.pow(6*self.spacing-math.pi*self.height/4,3)/self.length_skin - (0.5*math.pow(5*self.spacing-math.pi*self.height/4,2) - 1/6 * math.pow(5*self.spacing-math.pi*self.height/4,3)/self.length_skin)) + self.sum_boom_SC_y(5,5)*self.spacing) + self.shear_flow_magnitude_y[5]*self.spacing
		self.shear_flow_integrated_y[7] = -1/self.I_zz*	(self.skin_thickness*self.height/2*(0.5*math.pow(6.5*self.spacing-math.pi*self.height/4,2) - 1/6 * math.pow(6.5*self.spacing-math.pi*self.height/4,3)/self.length_skin - (0.5*math.pow(6*self.spacing-math.pi*self.height/4,2) - 1/6 * math.pow(6*self.spacing-math.pi*self.height/4,3)/self.length_skin)) + self.sum_boom_SC_y(6,6)*0.5*self.spacing) + self.shear_flow_magnitude_y[6]*0.5*self.spacing
		self.int3 = self.shear_flow_integrated_y[2] + self.shear_flow_integrated_y[3] + self.shear_flow_integrated_y[4] + self.shear_flow_integrated_y[5] + self.shear_flow_integrated_y[6] + self.shear_flow_integrated_y[7]
		if self.plane == "B737":
			self.shear_flow_integrated_y[7] = -1/self.I_zz*	(self.skin_thickness*self.height/2*(0.5*math.pow(7*self.spacing-math.pi*self.height/4,2) - 1/6 * math.pow(7*self.spacing-math.pi*self.height/4,3)/self.length_skin - (0.5*math.pow(6*self.spacing-math.pi*self.height/4,2) - 1/6 * math.pow(6*self.spacing-math.pi*self.height/4,3)/self.length_skin)) + self.sum_boom_SC_y(6,6)*self.spacing) + self.shear_flow_magnitude_y[6]*self.spacing
			self.shear_flow_integrated_y[8] = -1/self.I_zz*	(self.skin_thickness*self.height/2*(0.5*math.pow(7.5*self.spacing-math.pi*self.height/4,2) - 1/6 * math.pow(7.5*self.spacing-math.pi*self.height/4,3)/self.length_skin - (0.5*math.pow(7*self.spacing-math.pi*self.height/4,2) - 1/6 * math.pow(7*self.spacing-math.pi*self.height/4,3)/self.length_skin)) + self.sum_boom_SC_y(7,7)*0.5*self.spacing) + self.shear_flow_magnitude_y[7]*0.5*self.spacing
			self.int3 = self.shear_flow_integrated_y[2] + self.shear_flow_integrated_y[3] + self.shear_flow_integrated_y[4] + self.shear_flow_integrated_y[5] + self.shear_flow_integrated_y[6] + self.shear_flow_integrated_y[7] + self.shear_flow_integrated_y[8]

		#Region Y 4:
		self.shear_flow_integrated_y[8+self.plane_delta] = -1/self.I_zz*	(self.skin_thickness*self.height/2/self.length_skin*-1/6*math.pow(0.5*self.spacing,3)) + self.shear_flow_magnitude_y[7+self.plane_delta]*0.5*self.spacing
		self.shear_flow_integrated_y[9+self.plane_delta] = -1/self.I_zz*	(self.skin_thickness*self.height/2/self.length_skin*-1/6*(math.pow(1.5*self.spacing,3) - math.pow(0.5*self.spacing,3)) + self.sum_boom_SC_y(7+self.plane_delta,7+self.plane_delta))*self.spacing + self.shear_flow_magnitude_y[8+self.plane_delta]*self.spacing
		self.shear_flow_integrated_y[10+self.plane_delta] = -1/self.I_zz*	(self.skin_thickness*self.height/2/self.length_skin*-1/6*(math.pow(2.5*self.spacing,3) - math.pow(1.5*self.spacing,3)) + self.sum_boom_SC_y(8+self.plane_delta,8+self.plane_delta))*self.spacing + self.shear_flow_magnitude_y[9+self.plane_delta]*self.spacing
		self.shear_flow_integrated_y[11+self.plane_delta] = -1/self.I_zz*	(self.skin_thickness*self.height/2/self.length_skin*-1/6*(math.pow(3.5*self.spacing,3) - math.pow(2.5*self.spacing,3)) + self.sum_boom_SC_y(9+self.plane_delta,9+self.plane_delta))*self.spacing + self.shear_flow_magnitude_y[10+self.plane_delta]*self.spacing
		self.shear_flow_integrated_y[12+self.plane_delta] = -1/self.I_zz*	(self.skin_thickness*self.height/2/self.length_skin*-1/6*(math.pow(4.5*self.spacing,3) - math.pow(3.5*self.spacing,3)) + self.sum_boom_SC_y(10+self.plane_delta,10+self.plane_delta))*self.spacing + self.shear_flow_magnitude_y[11+self.plane_delta]*self.spacing
		self.shear_flow_integrated_y[13+self.plane_delta] = -1/self.I_zz*	(self.skin_thickness*self.height/2/self.length_skin*-1/6*(math.pow(self.length_skin,3) - math.pow(4.5*self.spacing,3)) + self.sum_boom_SC_y(11,11))*(self.length_skin-4.5*self.spacing) + self.shear_flow_magnitude_y[12]*(self.length_skin-4.5*self.spacing)
		self.int4 = self.shear_flow_integrated_y[8] + self.shear_flow_integrated_y[9] + self.shear_flow_integrated_y[10] + self.shear_flow_integrated_y[11] + self.shear_flow_integrated_y[12] + self.shear_flow_integrated_y[13]
		if self.plane == "B737":
			self.shear_flow_integrated_y[14] = -1/self.I_zz*	(self.skin_thickness*self.height/2/self.length_skin*-1/6*(math.pow(5.5*self.spacing,3) - math.pow(4.5*self.spacing,3)) + self.sum_boom_SC_y(12,12))*self.spacing + self.shear_flow_magnitude_y[13]*self.spacing
			self.shear_flow_integrated_y[15] = -1/self.I_zz*	(self.skin_thickness*self.height/2/self.length_skin*-1/6*(math.pow(self.length_skin,3) - math.pow(5.5*self.spacing,3)) + self.sum_boom_SC_y(13,13))*(self.length_skin-5.5*self.spacing) + self.shear_flow_magnitude_y[14]*(self.length_skin-5.5*self.spacing)
			self.int4 = self.shear_flow_integrated_y[9] + self.shear_flow_integrated_y[10] + self.shear_flow_integrated_y[11] + self.shear_flow_integrated_y[12] + self.shear_flow_integrated_y[13] + self.shear_flow_integrated_y[14] + self.shear_flow_integrated_y[15]
		#Region Y 5:
		self.shear_flow_integrated_y[16+2*self.plane_delta] = 1/self.I_zz *0.5*self.spar_thickness*math.pow(self.height/2,3)/3
		self.int5 = self.shear_flow_integrated_y[16]
			
		#Region Y 6:
		self.shear_flow_integrated_y[14+2*self.plane_delta] = -1/self.I_zz*	((-math.sin(-self.spacing/(self.height/2))+math.sin(-math.pi/2))*self.skin_thickness*math.pow(self.height/2,3)) + (self.shear_flow_magnitude_y[13+2*self.plane_delta]-self.shear_flow_magnitude_y[16+2*self.plane_delta])*(-self.spacing/(self.height/2)+math.pi/2)
		self.shear_flow_integrated_y[15+2*self.plane_delta] = -1/self.I_zz*	((-math.sin(0)+math.sin(-self.spacing/(self.height/2)))*self.skin_thickness*math.pow(self.height/2,3)+self.sum_boom_SC_y(12+2*self.plane_delta,12+2*self.plane_delta)*(self.spacing/(self.height/2)))+self.shear_flow_magnitude_y[14+2*self.plane_delta]*(self.spacing/(self.height/2))
		self.int6 = self.shear_flow_integrated_y[14+2*self.plane_delta] + self.shear_flow_integrated_y[15+2*self.plane_delta]
		return
			
	def qs0(self):#in this part of the code the only thing that needs to be added is the sum of the shear flows through the arc, sum of the shear flows through the straight part of the skin and the shear flow through the spar
		self.integrate_shear_flows()
		radius_arc = self.height/2 #defining the radius of the front section
		length_straight_skin = 2*self.length_skin 
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
								[-self.height/self.spar_thickness, +self.height/self.spar_thickness+2*self.length_skin/self.skin_thickness, -self.G*2*self.a2]])

		D = np.matrix([[1],[0],[0]])

		qs0 = np.linalg.solve(A,B)
		J_matrix = np.linalg.solve(C,D)
		self.J = 1/(J_matrix[2,0]*self.G)

		self.qs01 = float(qs0[0])
		self.qs02 = float(qs0[1])
		print("Q01:", self.qs01, "QS02:", self.qs02)
		print("Torsional stiffness:", self.J)
		return

	def shear_center(self):
		self.qs0()
		shear_center_z = (self.int1+self.int6)*self.height/2 + (self.int3+self.int4)* (self.height/2/self.length_skin) * (self.chord-self.height/2) + 2*self.a1*self.qs01 + 2*self.a2*self.qs02 - self.height/2
		print(self.height/2)
		return shear_center_z

	def normal_stress(self, z_pos, y_pos, Moment_z = 0, Moment_y = 0):
		if z_pos == 0:
			thickness = self.spar_thickness
		else:
			thickness = self.skin_thickness
		normal_stress = Moment_y*(z_pos-self.centroid_z)/self.I_yy	+ Moment_z*y_pos/self.I_zz
		return normal_stress
		
	def compute_shear_stress(self, y_shear, z_shear):
		self.actual_shear_flow()
		for i in range(0,len(self.shear_flow_magnitude_y)-2):
			self.shear_stress[i] = (self.shear_flow_magnitude_y[i]*y_shear + self.shear_flow_magnitude_z[i]*z_shear)/self.skin_thickness
		self.shear_stress[len(self.shear_flow_magnitude_y)-2] = (self.shear_flow_magnitude_y[len(self.shear_flow_magnitude_y)-2]*y_shear + self.shear_flow_magnitude_z[len(self.shear_flow_magnitude_y)-2]*z_shear)/self.skin_thickness
		self.shear_stress[len(self.shear_flow_magnitude_y)-1] = (self.shear_flow_magnitude_y[len(self.shear_flow_magnitude_y)-1]*y_shear + self.shear_flow_magnitude_z[len(self.shear_flow_magnitude_y)-1]*z_shear)/self.skin_thickness
		return self.shear_stress

	def Von_misses_stress(self, shear_stress, normal_stress):
		return math.sqrt(math.pow(normal_stress,2)+3*math.pow(shear_stress,2))
	
	def Compute_section(self, Moment_z, Moment_y, Shear_z, Shear_y):
		Moment_z = 5
		Moment_y = 5
		Shear_z = 5
		Shear_y	= 5
		
		z_circle = np.linspace(-self.height/2, 0, 1500)
		y_circle = np.sqrt((-self.height/2)**2 - (z_circle)**2)
		z_circle_2 = np.linspace(0, -self.height/2, 1500)
		y_circle_2 = -np.sqrt((-self.height/2)**2 - (z_circle_2)**2)
		z_plate = np.linspace(0, self.chord - self.height/2, 1500)
		y_plate = self.height/2 - (self.height/2)/(self.chord-self.height/2) * z_plate
		z_plate_2 = np.linspace(self.chord - self.height/2, 0, 1500)
		y_plate_2 = -self.height/2 + (self.height/2)/(self.chord-self.height/2) * z_plate_2
		y_vplate = np.linspace(-self.height/2, self.height/2,500)
		z_vplate = y_vplate*0

		z_pos = np.append(np.append(np.append(np.append(z_circle, z_plate), z_plate_2), z_circle_2),z_vplate)
		y_pos = np.append(np.append(np.append(np.append(y_circle, y_plate), y_plate_2), y_circle_2),y_vplate)
		
		von_misses = np.zeros(len(z_pos))

		normal_stress = np.zeros(len(z_pos))
		for i in range(0,len(normal_stress)):
			normal_stress[i] = self.normal_stress(z_pos[i], y_pos[i], Moment_z, Moment_y)
		self.compute_shear_stress(Shear_y, Shear_z)

		transversed = 0
		shear_index = 0
		for i in range(0,len(z_circle)):
			if z_circle[i] > self.booms_z[shear_index+1]:
				shear_index = shear_index + 1
				print(" Increasing +1: top circle",shear_index)
				
			von_misses[i] = self.Von_misses_stress(self.shear_stress[shear_index], normal_stress[i])

		shear_index = 2
		for i in range(0,len(z_plate)):
			if z_plate[i] > self.booms_z[shear_index] and self.booms_z[shear_index] > self.booms_z[shear_index-1]:
				shear_index = shear_index+1
				print(" Increasing +1: top plate",shear_index)
				
			von_misses[i+1500] = self.Von_misses_stress(self.shear_stress[shear_index], normal_stress[i])

		shear_index = shear_index + 1
		for i in range(0,len(z_plate_2)):
			if z_plate_2[i] < self.booms_z[shear_index-1]:
				shear_index = shear_index+1
				print(" Increasing bottom plate",shear_index)
				
			von_misses[i+3000] = self.Von_misses_stress(self.shear_stress[shear_index], normal_stress[i])

		shear_index = shear_index+1
		for i in range(0,len(z_circle_2)):
			if shear_index < (len(self.booms_z)+2) and z_circle_2[i] < self.booms_z[shear_index-2]:
				shear_index = shear_index+1
				print(" Increasing +1: bottom circle",shear_index)
			von_misses[i+4500] = self.Von_misses_stress(self.shear_stress[shear_index], normal_stress[i])

		shear_index = shear_index+1
		for i in range(0,len(z_vplate)):
			if y_vplate[i] > 0 and y_vplate[i-1] < 0:
				shear_index = shear_index+1
				print(" Increasing +1: vertical plate",shear_index)
				
			von_misses[i+6000] = self.Von_misses_stress(self.shear_stress[shear_index], normal_stress[i])
	

		von_mises_max = np.max(von_misses)
		von_mises_min = np.min(von_misses)
		fig = plt.figure()
		ax = fig.add_subplot(111)
		#x_boom, y_boom = self.booms_z, self.booms_y
		#ax.scatter(x_boom, y_boom)
		#ax.scatter(self.centroid_z, self.centroid_y)
		img = ax.scatter(z_pos, y_pos, c=von_misses, cmap=plt.cm.RdYlGn) #cmap = ((255*(von_misses[i]-von_mises_min)/(von_mises_max-von_mises_min)),(255*(1-(von_misses[i]-von_mises_min)/(von_mises_max-von_mises_min))),0))
		fig.colorbar(img)
		self.shear_center_z = self.shear_center()
		#ax.scatter(self.shear_center_z,0)
		ax.set_aspect(aspect=1)
		
		fig = plt.figure()
		ax = fig.add_subplot(111)
		ax.plot(z_pos,von_misses)
		plt.show()




x = Geometry(17.3/100,1.1/1000,2.5/1000,1.2/1000,1.4/100,1.8/100,0.484,13,1, "CRJ700")
#x.integrate_shear_flows()
#print("Shear center location CRJ700 is:*******************************************************", x.shear_center())
x.Compute_section(5, 5, 2, 2)

z =  [-0, -0.03692049, -0.12081074, -0.20884515, -0.29687956, -0.38491397, -0.47294838, -0.56098279, -0.56098279, -0.47294838, -0.38491397, -0.29687956, -0.20884515, -0.12081074 , -0.03692049]
z2 = []
z3 =  [-0.1025, -0.06557951, 0.018310740000000006, 0.10634515000000001, 0.19437956, 0.28241397, 0.37044838, 0.45848279, 0.45848279, 0.37044838, 0.28241397, 0.19437956, 0.10634515000000001, 0.018310740000000006, -0.06557951]


b737 = Geometry(20.5/100,1.1/1000,2.8/1000,1.2/1000,1.6/100,1.9/100,0.605,15,1, "B737")
#b737.idealization()
#print("Shear center location B737 is:*******************************************************", b737.shear_center())

#for a in z:
#	z2.append(-a-b737.height/2)
#print("THE LIST Z2 IS:", z2)
#print("Shear center location is:", x.shear_center())
