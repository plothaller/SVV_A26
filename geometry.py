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

'''


class Geometry:


	def __init__(self, height, skin_t, thickness_str, height_str, width_str, chord, number_str):
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
		self.booms_z = self.booms(self.spacing)[0]
		self.booms_y = self.booms(self.spacing)[1]
		self.zero_area_boom_z = [0,0,chord-height/2]
		self.zero_area_boom_y = [height/2,-height/2,0]
		self.SNx = self.booms(self.spacing/2)[0][1::2]
		self.SNy = self.booms(self.spacing/2)[1][1::2]
		self.I_zz = self.moments_of_inertia(self.booms_z, self.booms_y)[0]
		self.I_yy = self.moments_of_inertia(self.booms_z, self.booms_y)[1]

	def booms(self, spacing):
		print("Running booms")
		print("spacing:", spacing)
		n_str = self.spacing/spacing * self.n_str
		print("perimeter:", self.perimeter)
		x = []
		y = []
		if spacing > (math.pi*self.height/4):
			raise ValueError('The spacing is larger than the quarter circle')
		for i in range(0,math.ceil(n_str/2)):
			print(i)
			effective_spacing = spacing*i
			if effective_spacing == 0:
				print("=0", i)
				x.append(-self.height/2)
				y.append(0)
				incircle = 1
			elif effective_spacing < (math.pi*self.height/4):
				print("<pi*r/2", i)
				print("effective spacing is:", effective_spacing)
				incircle = incircle + 1
				x.append(-math.cos(2*effective_spacing/self.height)*self.height/2)
				y.append(math.sin(2*effective_spacing/self.height)*self.height/2)
			else:
				print(">pi*r/2", i)
				x.append(((effective_spacing-(math.pi*self.height/4))*((self.chord-self.height/2)/(math.sqrt(math.pow((self.chord-self.height/2),2)+math.pow((self.height/2),2))))))
				y.append(self.height/2*(1-((effective_spacing-(math.pi*self.height/4))*((self.chord-self.height/2)/(math.sqrt(math.pow((self.chord-self.height/2),2)+math.pow((self.height/2),2)))))/(self.chord-self.height/2)))
		for i in range(1, math.ceil(n_str/2)):
			x.append(x[i])
			y.append(-y[i])
		return x,y

	def idealization(self):
		fig = plt.figure()
		ax = fig.add_subplot(111)
		x_boom, y_boom = self.booms(self.spacing)
		x_circle = np.linspace(-self.height/2, 0, 50)
		y_circle = np.sqrt((-self.height/2)**2 - (x_circle)**2)
		x_plate = np.linspace(0, self.chord - self.height/2, 50)
		y_plate = self.height/2 - (self.height/2)/(self.chord-self.height/2) * x_plate
		ax.plot(x_circle, y_circle,'b')
		ax.plot(x_circle, -y_circle,'b')
		ax.plot(x_plate, y_plate,'b')
		ax.plot(x_plate, -y_plate,'b')
		ax.scatter(x_boom, y_boom)
		ax.scatter(self.SNx, self.SNy)
		ax.scatter(self.zero_area_boom_z, self.zero_area_boom_y)
		ax.set_aspect(aspect=1)
		plt.show()
	
	def centroid(self):
		print("No")

	def moments_of_inertia(self, z_boom, y_boom):
		I_zz = 0
		I_yy = 0
		for z in z_boom:
			I_yy = I_yy + math.pow(abs(z),2) * self.str_area
		for y in y_boom:
			I_zz = I_zz + math.pow(abs(y),2) * self.str_area 
		return I_zz, I_yy


	#def node_closest_to_x0(self, x_boom, y_boom):
	#	min_boom = min(x_boom, key=abs)
	#	min_boom_index = x_boom.index(min(x_boom, key=abs))
	#	print("Index is:", min_boom_index)

	def shear_center_x(self):
		shear_nodes_x = np.asarray(self.SNx)
		shear_nodes_y = np.asarray(self.SNy)
		shear_nodes_flow_z = np.zeros(len(self.booms_y))
		shear_nodes_flow_y = np.zeros(len(self.booms_y))
		Delta_theta = 0
		Delta_lenght = math.sqrt(math.pow(min(i for i in self.booms_z if i > 0),2)+math.pow(self.height/2 - self.booms_y[self.booms_z.index(min(i for i in self.booms_z if i > 0))],2))

		for i in range(0,len(self.SNx)):
			if -self.height/2 < self.booms_z[i] < 0: #Cases 1 and 6
				print("Doing shear center for node:", i, "Case 1,6")
				Delta_theta = Delta_theta + self.height/self.spacing
				if Delta_theta > math.pi/2:
					Delta_theta = math.pi/2
				shear_flow = -1/self.I_zz * ((math.cos(Delta_theta))*(-self.skin_thickness * math.pow(self.height,2))+(self.sum_booms_SC(0,i)))
				shear_nodes_flow_z[i] = 0
				print(shear_flow)
			else: 
				Delta_theta = 0
			if self.booms_z[i] > 0 and self.booms_y[i] > 0: #Case 3
				print("Doing shear center for node:", i, "Case 3")
				Delta_lenght = Delta_lenght + self.spacing
				shear_flow = -1/self.I_zz * (self.skin_thickness*self.height*Delta_lenght - (math.pow(Delta_lenght,2)/2*(self.skin_thickness*self.height)/((self.perimeter - math.pi*self.height/2)/2)))
				shear_nodes_flow_z[i] = 0
				print(shear_flow)
			elif self.booms_z[i] > 0 and self.booms_y[i] < 0: #Case 4
				print("Doing shear center for node:", i, "Case 4")
				Delta_lenght = Delta_lenght + self.spacing
				Shear_flow = -(math.pow(Delta_lenght,2)/2*(self.skin_thickness*self.height)/((self.perimeter - math.pi*self.height/2)/2))
				shear_nodes_flow_z[i] = 0
				print("ashas",shear_flow)
			else:
				Delta_lenght = math.sqrt(math.pow(min(i for i in self.booms_z if i > 0),2)+math.pow(self.height/2 - self.booms_y[self.booms_z.index(min(i for i in self.booms_z if i > 0))],2))

			
	def sum_booms_SC(self, start, end):
		summation = 0
		if end > len(self.booms_y):
			raise ValueError('Sum_boom_areas. End point is greater than the number of y_booms')
		for i in range(0, end):
			summation = summation + abs(self.booms_y[i]) * self.str_area
		return summation


x = Geometry(10,1,6,7,8,40,29)
print(x.idealization())
x_booms, y_booms = x.booms(x.spacing)
print("str_area:", x.str_area)
print("str 1:", y_booms[1])
print("Sum booms SC:", x.sum_booms_SC(0,2))
x.shear_center_x()

