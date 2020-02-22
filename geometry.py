import matplotlib.pyplot as plt
import numpy as np
import math


'''
TODO:

-Coordinate system uses X-Y (code) instead on Z-Y (report)
-Shear center
-Centroid
-Fix boom spacing
-What happens when the boom is exactly in the vertical plate and leading edge and trailing edge plate joint?????????
-Trailing edge shear element is not considered
-Aliron height is the full spar or only half the spar?????
-Shear flow cases 2,5 (spar) are yet not implemented
-Do we want to implement skin into the booms with a large number of booms (10k-100k) ????
-Shear flows are wrong, they are just the integral component
-Shear flow calculation uses wrong integration bounds

'''


class Geometry:


	def __init__(self, height, skin_t, thickness_str, height_str, width_str, chord, number_str, booms_per_str):
		self.height = height			#height of aileron
		self.skin_thickness = skin_t	#Thickness of skin
		self.t_st = thickness_str		#thickness of stringer
		self.h_st = height_str			#height of stringer
		self.w_st = width_str			#width of stringer
		self.chord = chord				#chord length
		self.n_str = number_str			#number of stringers
		self.str_area = width_str * thickness_str + height_str * thickness_str
		self.perimeter = math.pi*height/2 + 2*math.sqrt(math.pow(height/2,2) + math.pow(chord - height/2,2))
		self.spacing = (math.pi*height/2 + 2*math.sqrt(math.pow(height/2,2) + math.pow(chord - height/2,2)))/(number_str)
		self.booms_per_str = booms_per_str
		self.booms_z, self.booms_y, self.booms_area = self.booms(self.booms_per_str, True)[0], self.booms(self.booms_per_str, True)[1], self.booms(self.booms_per_str, True)[2]
		self.SNx = np.array(self.booms(self.booms_per_str*2)[0][1::2]) #np.append(self.booms(self.booms_per_str*2)[0][1::2], np.array([0, 0, (self.chord - self.height/2) - 1/4 * self.spacing, (self.chord - self.height/2) - 1/4 * self.spacing ]))
		self.SNy = np.array(self.booms(self.booms_per_str*2)[1][1::2]) #np.append(self.booms(self.booms_per_str*2)[1][1::2], np.array([self.height/4, -self.height/4 , self.height/2*(1-((self.chord - self.height/2) - 1/4 * self.spacing)/(self.chord-self.height/2)), -self.height/2*(1-((self.chord - self.height/2) - 1/4 * self.spacing)/(self.chord-self.height/2))]))
		self.centroid_z = self.centroid()[0]
		self.centroid_y = self.centroid()[1]
		self.I_zz = self.moments_of_inertia(self.booms_z, self.booms_y)[0]
		self.I_yy = self.moments_of_inertia(self.booms_z, self.booms_y)[1]

	def booms(self, booms_per_str, Booms = False):
		print("Running booms")
		spacing = self.spacing/booms_per_str
		n_str = self.n_str * booms_per_str
		vertical_plate_n_str = math.ceil(self.height/(2*spacing))
		print("spacing:", spacing)
		print("perimeter:", self.perimeter)
		x = []
		y = []
		a = []
		if spacing > (math.pi*self.height/4):
			raise ValueError('The spacing is larger than the quarter circle')
		for i in range(0,math.ceil(n_str/2)):
			print(i)
			effective_spacing = spacing*i
			if effective_spacing == 0: #Most fordward point (L.E line)
				print("=0", i)
				x.append(-self.height/2)
				y.append(0)
				if i%booms_per_str == 0:
					a.append(self.str_area)
				else:
					a.append(0)
			elif effective_spacing < (math.pi*self.height/4): #Points in the circle (except the most front one)
				print("<pi*r/2", i)
				print("effective spacing is:", effective_spacing)
				x.append(-math.cos(2*effective_spacing/self.height)*self.height/2)
				y.append(math.sin(2*effective_spacing/self.height)*self.height/2)
				if i%booms_per_str == 0:
					a.append(self.str_area)
				else:
					a.append(0)
			else:	#Points in the straight part 
				print(">pi*r/2", i)
				x.append(((effective_spacing-(math.pi*self.height/4))*((self.chord-self.height/2)/(math.sqrt(math.pow((self.chord-self.height/2),2)+math.pow((self.height/2),2))))))
				y.append(self.height/2*(1-((effective_spacing-(math.pi*self.height/4))*((self.chord-self.height/2)/(math.sqrt(math.pow((self.chord-self.height/2),2)+math.pow((self.height/2),2)))))/(self.chord-self.height/2)))
				if i%booms_per_str == 0:
					a.append(self.str_area)
				else:
					a.append(0)
		for i in range(1, len(x)):
			x.append(x[i])
			y.append(-y[i])
			a.append(a[i])
		if Booms == True:
			 x.append(0)
			 y.append(-self.height/2)
			 a.append(0)
		for i in range(1,(vertical_plate_n_str)):
			x.append(0)
			y.append(-self.height/2 +i/vertical_plate_n_str*self.height)
			a.append(0)
		if Booms == True:
			 x.append(0)
			 y.append(self.height/2)
			 a.append(0)

		return x,y,a

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
		ax.scatter(self.SNx, self.SNy)
		ax.scatter(self.centroid_z, self.centroid_y)
		ax.set_aspect(aspect=1)
		plt.show()
	
	def centroid(self):
		centroid_z = 0
		centroid_y = 0
		for z in self.booms_z:
			centroid_z = centroid_z + z 
		return centroid_z/len(self.booms_z), centroid_y

	def moments_of_inertia(self, z_boom, y_boom):
		I_zz = 0
		I_yy = 0
		for z in z_boom:
			I_yy = I_yy + math.pow(abs(z-self.centroid_z),2) * self.str_area
		for y in y_boom:
			I_zz = I_zz + math.pow(abs(y-self.centroid_y),2) * self.str_area 
		return I_zz, I_yy

	def shear_center_x(self):
		print("np.zeros(1):", np.zeros(1))
		shear_nodes_x = np.asarray(self.SNx)
		shear_nodes_y = np.asarray(self.SNy)
		shear_nodes_flow_z = np.zeros(self.SNx.size)
		shear_nodes_flow_y = np.zeros(self.SNx.size)
		Delta_theta = 0
		Delta_lenght = math.sqrt(math.pow(min(i for i in self.booms_z if i > 0),2)+math.pow(self.height/2 - self.booms_y[self.booms_z.index(min(i for i in self.booms_z if i > 0))],2))

		for i in range(0,len(self.SNx)):
			if -self.height/2 < self.SNx[i] < 0: #Cases 1 and 6
				print("Doing shear center for node:", i, "Case 1,6")
				Delta_theta = Delta_theta + self.height/self.spacing
				if Delta_theta > math.pi/2:
					Delta_theta = math.pi/2
				shear_flow = -1/self.I_zz * ((math.cos(Delta_theta))*(-self.skin_thickness * math.pow(self.height,2))+(self.sum_booms_SC(0,i)))
				shear_nodes_flow_z[i] = 0
				print(shear_flow)
			else: 
				Delta_theta = 0
			if self.SNx[i] > 0 and self.SNy[i] > 0: #Case 3
				print("Doing shear center for node:", i, "Case 3")
				Delta_lenght = Delta_lenght + self.spacing
				#The shear flows are wrong, they are only the integral component
				shear_flow = -1/self.I_zz * (self.skin_thickness*self.height*Delta_lenght - (math.pow(Delta_lenght,2)/2*(self.skin_thickness*self.height)/((self.perimeter - math.pi*self.height/2)/2)))
				shear_nodes_flow_z[i] = 0
				print(shear_flow)
			elif self.SNx[i] > 0 and self.SNy[i] < 0: #Case 4
				print("Doing shear center for node:", i, "Case 4")
				Delta_lenght = Delta_lenght + self.spacing
				#The shear flows are wrong, they are only the integral component
				Shear_flow = -(math.pow(Delta_lenght,2)/2*(self.skin_thickness*self.height)/((self.perimeter - math.pi*self.height/2)/2))
				shear_nodes_flow_z[i] = 0
				print("ashas",shear_flow)
			else:
				Delta_lenght = math.sqrt(math.pow(min(i for i in self.booms_z if i > 0),2)+math.pow(self.height/2 - self.booms_y[self.booms_z.index(min(i for i in self.booms_z if i > 0))],2))
			if self.SNx[i] == 0 and self.SNy[i] == 0:
				Delta_lenght = self.height/2
				#The shear flows are wrong, they are only the integral component
	
	def sum_booms_SC(self, start, end):
		summation = 0
		if end > len(self.booms_y):
			raise ValueError('Sum_boom_areas. End point is greater than the number of y_booms')
		for i in range(0, end):
			summation = summation + abs(self.booms_y[i]) * self.str_area
		return summation


x = Geometry(10,1,6,7,8,40,29,1)
print(x.idealization())
x_booms, y_booms = x.booms(x.spacing)[0], x.booms(x.spacing)[1] 
print("str_area:", x.str_area)
print("str 1:", y_booms[1])
print("Sum booms SC:", x.sum_booms_SC(0,2))
print(x.shear_center_x())

