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
        x_circle = np.linspace(-self.h/2, 0, 50)
        y_circle = np.sqrt((-self.h/2)**2 - (x_circle)**2)
        y_circle = np.append(y_circle, -y_circle)
        plt.plot(x_circle, y_circle,'b')
        plt.plot(x_circle, -y_circle)
        plt.show()


    def area_str(self):
        self.area_st = (self.h_st + self.w_st)* self.t_st
        return self.area_st


x = Geometry(5,6,7,8,9,10)
print(x.idealization())

