import matplotlib.pyplot as plt
import numpy as np
import math


'''
TODO:

-Coordinate system uses X-Y (code) instead on Z-Y (report)
-Shear center
-Centroid
-Fix boom spacing

'''


class Geometry:

	def __init__(self, height, thickness_str, height_str, width_str, chord, number_str):
		self.height = height        #height of aileron
		self.t_st = thickness_str   #thickness of stringer
		self.h_st = height_str      #height of stringer
		self.w_st = width_str       #width of stringer
		self.chord = chord          #chord length
		self.n_str = number_str     #number of stringers
		self.str_area = width_str * thickness_str + height_str * thickness_str
		self.perimeter = math.pi*height/2 + 2*math.sqrt(math.pow(height/2,2) + math.pow(chord - height/2,2))
		self.spacing = (math.pi*height/2 + 2*math.sqrt(math.pow(height/2,2) + math.pow(chord - height/2,2)))/(number_str)
		self.Stratos_scaling_factor = (self.height/2)/(self.chord-self.height/2) 

	def idealization(self):
		fig = plt.figure()
		ax = fig.add_subplot(111)
		x_boom, y_boom = self.booms()
		x_circle = np.linspace(-self.height/2, 0, 50)
		y_circle = np.sqrt((-self.height/2)**2 - (x_circle)**2)
		x_plate = np.linspace(0, self.chord - self.height/2, 50)
		y_plate = self.height/2 - (self.height/2)/(self.chord-self.height/2) * x_plate
		ax.plot(x_circle, y_circle,'b')
		ax.plot(x_circle, -y_circle,'b')
		ax.plot(x_plate, y_plate,'b')
		ax.plot(x_plate, -y_plate,'b')
		ax.scatter(x_boom, y_boom)
		ax.set_aspect(aspect=1)
		plt.show()
	
	def booms(self):
		print("Running booms")
		print("spacing:", self.spacing)
		print("perimeter:", self.perimeter)
		x = []
		y = []
		if self.spacing > (math.pi*self.height/4):
			raise ValueError('The spacing is larger than the quarter circle')
		for i in range(0,math.ceil(self.n_str/2)):
			print(i)
			effective_spacing = self.spacing*i
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
		for i in range(1, math.ceil(self.n_str/2)):
			x.append(x[i])
			y.append(-y[i])
		return x,y




	def compute_geometric_properties(self):
		x_boom, y_boom = self.booms()
		self.moments_of_inertia(x_boom, y_boom)

	def centroid(self):
		print("No")

	def moments_of_inertia(self, x_boom, y_boom):
		I_xx = 0
		I_yy = 0
		for x in x_boom:
			I_yy = I_yy + math.pow(abs(x),2) * self.str_area
		for y in y_boom:
			I_xx = I_xx + math.pow(abs(y),2) * self.str_area 
		return I_xx, I_yy


	def node_closet_to_x0(self, x_boom, y_boom):
		min_boom = min(x_boom, key=abs)
		min_boom_index = x_boom.index(min(x_boom, key=abs))
		print("Index is:", min_boom_index)



		


x = Geometry(10,6,7,8,40,29)
print(x.idealization())
print(x.booms())
print(x.compute_geometric_properties())
x_booms, y_booms = x.booms()
print(x.node_closet_to_x0(x_booms, y_booms))

