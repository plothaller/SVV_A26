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
-Aileron height is the full spar or only half the spar?????
-Shear flow cases 2,5 (spar) are yet not implemented
-Do we want to implement skin into the booms with a large number of booms (10k-100k) ????
-Shear flows are wrong, they are just the integral component
-Shear flow calculation uses wrong integration bounds
-On the circle (L.E) to plates interface, take into accunt the lenght of the plate between the last boom in the circle and the vertical plate.
-Be careful with the inputs, the height is the full aliron thickness, not just half

'''


class Geometry:
	def __init__(self, height, skin_t, spar_t, thickness_str, height_str, width_str, chord, number_str, booms_per_str):
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
		self.spacing_extra_booms = self.spacing/booms_per_str
		self.booms_per_str = booms_per_str
		self.booms_z, self.booms_y, self.booms_area = self.booms(self.booms_per_str, True)[0], self.booms(self.booms_per_str, True)[1], self.booms(self.booms_per_str, True)[2]
		self.SNx = np.array(self.booms(self.booms_per_str*2)[0][1::2]) #np.append(self.booms(self.booms_per_str*2)[0][1::2], np.array([0, 0, (self.chord - self.height/2) - 1/4 * self.spacing, (self.chord - self.height/2) - 1/4 * self.spacing ]))
		self.SNy = np.array(self.booms(self.booms_per_str*2)[1][1::2]) #np.append(self.booms(self.booms_per_str*2)[1][1::2], np.array([self.height/4, -self.height/4 , self.height/2*(1-((self.chord - self.height/2) - 1/4 * self.spacing)/(self.chord-self.height/2)), -self.height/2*(1-((self.chord - self.height/2) - 1/4 * self.spacing)/(self.chord-self.height/2))]))
		self.centroid_z = self.centroid()[0]
		self.centroid_y = self.centroid()[1]
		self.I_zz = self.moments_of_inertia(self.booms_z, self.booms_y)[0]
		self.I_yy = self.moments_of_inertia(self.booms_z, self.booms_y)[1]
		self.Top_half = 0		#Assign values in shear center computation
		self.Top_plate = 0		#Assign values in shear center computation
		self.Bottom_half = 0	#Assign values in shear center computation
		self.Bottom_plate = 0	#Assign values in shear center computation
		self.shear_center_x = self.shear_center()[0]
		print("The booms are located in:", self.booms_z, self.booms_y)


	def idealization(self):
		x_circle = np.linspace(-self.h/2, 0, 500)
		y_circle = np.sqrt((-self.h/2)**2 - (x_circle)**2)
		x_skin = np.linspace(0, (self.c_a -self.h/2), 500)
		gradient_skin = ((-self.h/2)/(self.c_a - self.h/2))
		y_skin = gradient_skin*(x_skin) + self.h/2
		plt.plot(x_circle, y_circle,'b')
		plt.plot(x_skin, y_skin, 'b')
		plt.show()

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
		
		if Booms == False:
			x.append(-self.height/2)
			y.append(0)
			a.append(0)
		x.append(self.chord-self.height/2)
		y.append(0)
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
		#ax.scatter(self.shear_center_x,0)
		ax.set_aspect(aspect=1)
		print("The booms are located in		X: ", self.booms_z, "	Y:",self.booms_y)
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

	def shear_center(self):
		print("np.zeros(1):", np.zeros(1))
		shear_nodes_z = np.asarray(self.SNx)
		shear_nodes_y = np.asarray(self.SNy)
		shear_flow_magnitude = np.zeros(self.SNx.size)
		shear_nodes_flow_z = np.zeros(self.SNx.size)
		shear_nodes_flow_y = np.zeros(self.SNx.size)
		Delta_theta = 0
		Delta_lenght = math.sqrt(math.pow(min(i for i in self.booms_z if i > 0),2)+math.pow(self.height/2 - self.booms_y[self.booms_z.index(min(i for i in self.booms_z if i > 0))],2))
		self.Top_half = 0
		self.Top_plate = 0
		self.Bottom_half = 0
		self.Bottom_plate = 0
		skin_lenght_weighted_ratio = self.height/(4*math.sqrt(math.pow(self.chord-self.height/2,2)+math.pow(self.height/2,2)))
		x_to_s_region3 = (math.sqrt(math.pow(self.chord-self.height/2,2)+math.pow(self.height/2,2)))/(self.chord-self.height/2)


		#First we compute the Delta shear flow in each section
		for i in range(0,len(self.SNx)):
			print("X: ",self.SNx, "Y: ",self.SNy)
			if self.SNy[i] == 0:
				shear_flow_magnitude[i] = 0
				if self.Top_half !=0 and self.Top_plate == 0:
					self.Top_plate = i
				print("*******************************************************************************DOING SHEAR FLOW = 0 for i:", i)
			if -self.height/2 < self.SNx[i] < 0: #Cases 1 and 6
				print("Doing open shear center for node:", i, "Case 1,6")
				Delta_theta = Delta_theta + self.height/self.spacing
				shear_flow_magnitude[i] = -1/self.I_zz * (self.skin_thickness*math.pow(self.height/2,2)*(self.booms_z[i]-self.booms_z[i+1])/(self.height/2))
				Delta_theta = Delta_theta + self.spacing/(self.height/2)

			else:
				Delta_theta = 0
			if self.SNx[i] > 0 and self.SNy[i] > 0: #Case 3
				print("Doing open shear center for node:", i, "Case 3")
				if self.Top_half == 0:
					self.Top_half = i
				shear_flow_magnitude[i] = -1/self.I_zz *((self.height/2*self.skin_thickness)*(x_to_s_region3*self.booms_z[i+1]-x_to_s_region3*self.booms_z[i])+skin_lenght_weighted_ratio*self.skin_thickness*(math.pow(x_to_s_region3*self.booms_z[i],2)-math.pow(x_to_s_region3*self.booms_z[i+1],2)))
				Delta_lenght = Delta_lenght + self.spacing
			elif self.SNx[i] > 0 and self.SNy[i] < 0: #Case 4
				if self.Bottom_half == 0:
					self.Bottom_half = i
				print("Doing open shear center for node:", i, "Case 4")
				Delta_lenght = Delta_lenght + self.spacing
				#The shear flows are wrong, they are only the integral component
				shear_flow_magnitude[i] = 1/self.I_zz *((self.height/2*self.skin_thickness)*(x_to_s_region3*self.booms_z[i]-x_to_s_region3*self.booms_z[i+1])+skin_lenght_weighted_ratio*self.skin_thickness*(math.pow(x_to_s_region3*self.booms_z[i+1],2)-math.pow(x_to_s_region3*self.booms_z[i],2)))
			else:
				Delta_lenght = math.sqrt(math.pow(min(i for i in self.booms_z if i > 0),2)+math.pow(self.height/2 - self.booms_y[self.booms_z.index(min(i for i in self.booms_z if i > 0))],2))
			if self.SNx[i] == 0:
				print("Doing open shear center for node:", i, "Case 5,2")
				if self.Bottom_plate == 0:
					self.Bottom_plate = i
				shear_flow_magnitude[i] = -1/self.I_zz * (1/2*self.skin_thickness*(math.pow(self.booms_y[i],2)-math.pow(self.booms_y[i-1],2)))


				#The shear flows are wrong, they are only the integral component
		#self.booms_z = np.delete(self.booms_z, self.Top_plate)
		#self.booms_y = np.delete(self.booms_y, self.Top_plate)



		#self.SNx = np.delete(self.SNx, self.Bottom_plate)
		#self.SNy = np.delete(self.SNy, self.Bottom_plate)
		#shear_flow_magnitude = np.delete(shear_flow_magnitude, self.Bottom_plate)
		#shear_nodes_flow_y = np.delete(shear_nodes_flow_y, self.Bottom_plate)
		#shear_nodes_flow_z = np.delete(shear_nodes_flow_z, self.Bottom_plate)

		#Add previous node shear flow. Change some qs01 for qs02 !!!!!!!!!!!!!!!!!!!
		#Add the previus node shear flow!!!!!!!!!!!!!!!!!!!!!!
		qs01, qs02 = self.qs0(shear_flow_magnitude)
		for i in range(0,len(self.SNx)):
			print("Adding the constant closed cell shear flows q01 and q02 for node:", i, "	QS01:",qs01, "QS02:", qs02, "shear_flow_magnitude:", shear_flow_magnitude[i], "En eje z:", self.components(i)[0]*shear_flow_magnitude[i])
			if i < self.Top_half:
				shear_flow_magnitude[i] = shear_flow_magnitude[i] #- qs01
				print("AFTER DOING SMTH IT IS:",shear_flow_magnitude[i] )
			elif i < self.Top_plate:
				shear_flow_magnitude[i] = shear_flow_magnitude[i] #- qs02
			elif i < self.Bottom_half:
				self.booms_z[self.Top_plate+1] = self.SNx[self.Top_plate]
				self.booms_y[self.Top_plate+1] = self.SNy[self.Top_plate]
				print("HACIENDO LA PRIMERA DEL CIRC ABAJO:", i)
				shear_flow_magnitude[i] = shear_flow_magnitude[i] #+ qs01 #- shear_flow_magnitude[self.Top_half-1]
			elif i < self.Bottom_plate:
				shear_flow_magnitude[i] = shear_flow_magnitude[i] #+ qs02
			else:
				print("BOT: ",shear_flow_magnitude[self.Bottom_plate-3], "	TOP:",shear_flow_magnitude[self.Top_plate-1])
				difference =  shear_flow_magnitude[self.Bottom_plate-3] + shear_flow_magnitude[self.Top_plate-1]
				print("DIF:", difference)
				shear_flow_magnitude[i] = shear_flow_magnitude[i] + difference + shear_flow_magnitude[self.Bottom_half-1] #- qs02 + qs01	
				shear_nodes_flow_z[i] = 0
				shear_nodes_flow_y[i] = shear_flow_magnitude[i]
			#This next condition is wrong. SF at Top_half != 0. Apply the correct boundary condition!!!!!!!!!!!!!!!!!!!!!!!!  DO the shear in the vertical plate is equal to the sum of the shears and then the shear in the plate is the other stuff
			if i == self.Top_plate:
				shear_flow_magnitude[i] = 0
			elif i != 0 and i != self.Top_plate and i != self.Bottom_half and i != self.Bottom_plate:
				shear_flow_magnitude[i] = shear_flow_magnitude[i] -1/self.I_zz *shear_flow_magnitude[i-1]
			if abs(shear_flow_magnitude[i]) > 1:
				print("Shear magnitude too big:", i)
			if self.SNy[i] == 0:
				shear_flow_magnitude[i] = 0
			if shear_nodes_flow_y[i] == 0:
				shear_nodes_flow_z[i] = self.components(i)[0]*shear_flow_magnitude[i]
				shear_nodes_flow_y[i] = self.components(i)[1]*shear_flow_magnitude[i]




		fig = plt.figure()
		ax = fig.add_subplot(111)
		x_boom, y_boom = self.booms_z, self.booms_y
		x_circle = np.linspace(-self.height/2, 0, 50)
		y_circle = np.sqrt((-self.height/2)**2 - (x_circle)**2)
		x_plate = np.linspace(0, self.chord - self.height/2, 50)
		y_plate = self.height/2 - (self.height/2)/(self.chord-self.height/2) * x_plate
		y_vplate = np.linspace(-self.height/2, self.height/2, 2)
		x_vplate = y_vplate*0
		#ax.plot(x_circle, y_circle,'b')
		#ax.plot(x_circle, -y_circle,'b')
		#ax.plot(x_plate, y_plate,'b')
		#ax.plot(x_plate, -y_plate,'b')
		#ax.plot(x_vplate, y_vplate,'b')
		#ax.scatter(x_boom, y_boom)
		ax.scatter(self.SNx, self.SNy)
		ax.scatter(self.SNx[0], self.SNy[0])
		ax.scatter(self.SNx[self.Top_plate], self.SNy[self.Top_plate])
		#ax.scatter(self.SNx[self.Top_half], self.SNy[self.Top_half])
		#ax.scatter(self.SNx[self.Bottom_half], self.SNy[self.Bottom_half])
		ax.scatter(self.booms_z[self.Top_plate], self.booms_y[self.Top_plate])
		#ax.scatter(self.booms_z[self.Bottom_plate-1], self.booms_y[self.Bottom_plate-1])
		ax.scatter(self.booms_z[self.Top_plate+1], self.booms_y[self.Top_plate+1])

		ax.scatter(self.centroid_z, self.centroid_y)

		for i in range(0,len(self.SNx)):
			z_vect = [self.SNx[i], self.SNx[i] + shear_nodes_flow_z[i]/100000000000]
			y_vect = [self.SNy[i], self.SNy[i] + shear_nodes_flow_y[i]/100000000000]
			print("For index: ", i,"We have: ", z_vect, y_vect, "SF magnitude:", shear_flow_magnitude[i])
			ax.plot(z_vect, y_vect,'r')
		print("Boundaries:")
		print("Top_half:", self.Top_half)
		print("Top_plate:", self.Top_plate)
		print("Bottom_half:", self.Bottom_half)
		print("Bottom_plate:", self.Bottom_plate)
		print("Lenght:", len(self.SNx))
		#print("Shear_flow_magnitude i=3", shear_flow_magnitude[3])
		#print("Shear_flow_magnitude i=4", shear_flow_magnitude[4])
		#print("Shear_flow_magnitude i=17", shear_flow_magnitude[17])
		#print("Shear_flow_magnitude i=18", shear_flow_magnitude[18])
		#print("Shear_flow_magnitude i=29", shear_flow_magnitude[29])
		#print("Shear_flow_magnitude i=30", shear_flow_magnitude[30])


		ax.set_aspect(aspect=1)
		plt.show()


		#Multiply by the skin lenght
		shear_center_z = np.dot(shear_nodes_flow_z, shear_nodes_y) - np.dot(shear_nodes_flow_y, shear_nodes_z)
		print("The shear center is located at z:", shear_center_z)
		self.booms_z, self.booms_y, self.booms_area = self.booms(self.booms_per_str, True)[0], self.booms(self.booms_per_str, True)[1], self.booms(self.booms_per_str, True)[2]
		return shear_center_z, 0


	def components(self, index):
		Deltaz = self.booms_z[index+1] - self.booms_z[index]
		Deltay = self.booms_y[index+1] - self.booms_y[index]
		Lenght = math.sqrt(math.pow(Deltaz,2)+math.pow(Deltay,2))
		z_comp = Deltaz/Lenght
		y_comp = Deltay/Lenght
		return z_comp, y_comp
	
	#def sum_booms_SC(self, start, end):
	#	summation = 0
	#	if end > len(self.booms_y):
	#		raise ValueError('Sum_boom_areas. End point is greater than the number of y_booms')
	#	for i in range(0, end):
	#		summation = summation + self.booms_y[i] * self.str_area
	#	return summation


	def qs0(self, shear_flow_magnitude):#in this part of the code the only thing that needs to be added is the sum of the shear flows through the arc, sum of the shear flows through the straight part of the skin and the shear flow through the spar
		radius_arc = self.height/2 #defining the radius of the front section
		length_straight_skin = 2*np.sqrt(math.pow(radius_arc,2)+math.pow(self.chord - radius_arc,2))
		effective_lenght = 0
		sum_shearflow_through_arc = 0
		sum_shearflow_through_straightskin = 0
		shearflow_spar = 0
		for i in range(0,self.Top_half-1):
			effective_spacing = self.spacing_extra_booms
			if (i+1)*self.spacing_extra_booms > math.pi*self.height/4:
				effective_spacing = math.pi*self.height/4 - self.spacing_extra_booms
			if effective_spacing < 0:
				raise ValueError('effective lenght < 0 (arc, qs0)')
			sum_shearflow_through_arc = sum_shearflow_through_arc + shear_flow_magnitude[i] * effective_spacing
		for i in range(self.Top_plate,self.Bottom_half-1):
			effective_spacing = self.spacing_extra_booms
			if (i+1)*self.spacing_extra_booms > math.pi*self.height/4:
				effective_spacing = math.pi*self.height/4 - self.spacing_extra_booms
			if effective_spacing < 0:
				raise ValueError('effective lenght < 0 (arc, qs0)')
			sum_shearflow_through_arc = sum_shearflow_through_arc + shear_flow_magnitude[i] * effective_spacing 
		#sum_shearflow_through_arc = #multiply the shearflows with the distance they act over and sum them all for the arc part of the aileron (difficult if you get close to spar). Read the book, sent from my htc
		for i in range(self.Top_half, self.Top_plate-1):
			effective_lenght = self.spacing_extra_booms
			if self.SNx[i-1] < 0:
				effective_lenght = i*self.spacing_extra_booms - math.pi*self.height/4
			if effective_spacing < 0:
				raise ValueError('effective lenght < 0 (arc, qs0)')
			sum_shearflow_through_straightskin = sum_shearflow_through_straightskin + shear_flow_magnitude[i] * effective_spacing
		for i in range(self.Bottom_half, self.Bottom_plate-1):
			effective_lenght = self.spacing_extra_booms
			if self.SNx[i-1] < 0:
				effective_lenght = i*self.spacing_extra_booms - math.pi*self.height/4
			if effective_spacing < 0:
				raise ValueError('effective lenght < 0 (arc, qs0)')
			sum_shearflow_through_straightskin = sum_shearflow_through_straightskin + shear_flow_magnitude[i] * effective_spacing
		for i in range(self.Top_half, len(self.SNx)-1):
			effective_spacing = self.booms_y[i]-self.booms_y[i-1]
			shearflow_spar = sum_shearflow_through_straightskin + shear_flow_magnitude[i] * effective_spacing 
		shearflow_spar = 0#to be filled in by Guille. Multiplied with distance
		A = np.matrix([[np.pi*self.height/self.skin_thickness + self.height/self.spar_thickness, - self.height/self.spar_thickness],
                               [- self.height/self.spar_thickness, length_straight_skin/self.skin_thickness + self.height/self.spar_thickness]])
		B = np.matrix([[sum_shearflow_through_arc/self.skin_thickness - shearflow_spar/self.spar_thickness],
                               [sum_shearflow_through_straightskin/self.skin_thickness + shearflow_spar/self.spar_thickness]])
		qs0 = np.linalg.solve(A,B)
		qs01 = qs0[0]
		qs02 = qs0[1]
		return qs01, qs01

	#def qs0(self, shear_flow_magnitude):
	#	radius_arc = self.height/2 #defining the radius of the front section
	#	perp_dist_to_straight = (self.height/2 * (self.chord - (self.height / 2)))/(math.sqrt((self.height/2)**2) + (self.chord - (self.height / 2))**2) #perpendicular distance to the straight part of the aileron
	#	for i in range(0,len(self.SNx)):
	#		if -self.height/2 < self.SNx[i] < 0: #Cases 1 and 6
	#			qs01 = -(self.spacing * radius_arc * (np.sum(shear_flow_magnitude[0:self.Top_half-1])+np.sum(shear_flow_magnitude[self.Top_plate-1:self.Top_half-1]))) / (math.pi * math.pow(radius_arc, 2))
	#		if 0 < self.SNx[i] < (self.chord - self.height/2):
	#			qs02 = -(self.spacing * perp_dist_to_straight * (np.sum(shear_flow_magnitude[self.Top_half-1:self.Top_plate-1])+np.sum(shear_flow_magnitude[self.Bottom_half-1:self.Bottom_plate-1]))) / (2*radius_arc * (self.chord - radius_arc))
	#	return qs01, qs02


x = Geometry(17.3/100,1.1/1000,2.5/1000,1.2/1000,1.4/100,1.4/100,0.484,13,1)
x.idealization()
x_booms, y_booms = x.booms(x.spacing)[0], x.booms(x.spacing)[1] 
print("str_area:", x.str_area)
print("str 1:", y_booms[1])
print("Shear center location is:", x.shear_center())

