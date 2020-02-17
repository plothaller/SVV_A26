import matplotlib.pyplot as plt
import numpy as np


class Geometry:

	def __init__(self, height, thickness_str, height_str, width_str, chord, number_str):
		self.h = height             #height of aileron
		self.t_st = thickness_str   #thickness of stringer
		self.h_st = height_str      #height of stringer
		self.w_st = width_str       #width of stringer
		self.c_a = chord            #chord length
		self.n_st = number_str      #number of stringers

	def idealization(self):
		x_circle = np.linspace(-self.h/2, 0, 500)
		y_circle = np.sqrt((-self.h/2)**2 - (x_circle)**2)
		x_skin = np.linspace(0, (self.c_a -self.h/2), 500)
		gradient_skin = ((-self.h/2)/(self.c_a - self.h/2))
		y_skin = gradient_skin*(x_skin) + self.h/2
		plt.plot(x_circle, y_circle,'b')
		plt.plot(x_skin, y_skin, 'b')
		plt.show()


	def area_str(self):
		self.area_st = (self.h_st + self.w_st)* self.t_st
		return self.area_st


x = Geometry(0.173,0.0012,0.0014,0.0018,0.484,13)
print(x.idealization())

