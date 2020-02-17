import matplotlib.pyplot as plt
import numpy as np
import math


class Geometry:

	def __init__(self, height, thickness_str, height_str, width_str, chord, number_str):
		self.h = height             #height of aileron
		self.t_st = thickness_str   #thickness of stringer
		self.h_st = height_str      #height of stringer
		self.w_st = width_str       #width of stringer
		self.c_a = chord            #chord length
		self.n_str = number_str      #number of stringers
		self.perimeter = math.pi*height/2 + 2*math.sqrt(math.pow(height/2,2) + math.pow(chord - height/2,2))
		self.spacing = math.pi*height/2 + 2*math.sqrt(math.pow(height/2,2) + math.pow(chord - height/2,2))/number_str
		self.Stratos_scaling_factor = (height)/(2*chord) - 1 

	def idealization(self):
		x1_circle = np.linspace(-self.h/2, 0, 50)
		y1_circle = np.sqrt((-self.h/2)**2 - (x1_circle)**2)
		x_circle = np.append(x1_circle, x1_circle)
		y_circle = np.append(y1_circle, -y1_circle)
		#plt.plot(x_circle, y_circle,'b')
		plt.plot(x_circle, -y_circle)
		plt.show()


	def area_str(self):
		self.area_st = (self.h_st + self.w_st)* self.t_st
		return self.area_st
	
	def booms(self):
		x = []
		y = []
		if self.spacing > (math.pi*self.height/4):
			raise ValueError('The spacing is larger than the quarter circle')
		for i in range(0,self.n_str):
			if i*self.spacing <	(math.pi*self.height/4):
				x.append(self.height/2 * math.cos(self.spacing/(math.pi*self.height/2)*180))
				y.append(self.height/2 * math.sin(self.spacing/(math.pi*self.height/2)*180))
			else:
				x.append = math.sqrt((self.spacing-(math.pi*self.height)/4)/math.sqrt(math.pow(self.Stratos_scaling_factor,2)+1))
				y.append = self.Stratos_scaling_factor * x
		print(x,y)






x = Geometry(5,6,7,8,9,10)
print(x.idealization())

