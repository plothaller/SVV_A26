# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 14:01:53 2020

@author: Max v. Huffelen
"""

import numpy as np
import os

#global variables C_a and l_a for the chord and length of the aileron, respectively
C_a = 0.484
l_a = 1.691

def locationz(idx):
    #inputs: the index of the element in the matrix
    #returns the z-location of the element on the wing according to the equation given in the assignment
    idx += 1    #add the += 1 as python indexing begins at 0 rather than at 1
    #C_a = 1 #returns as a fraction of chord length. use the lenght of C_a in meters to get a value as such
    theta_z = lambda i : (i-1)*np.pi/81    
    z_prime = -0.5*(C_a/2*(1 - np.cos(theta_z(idx))) + C_a/2*(1 - np.cos(theta_z(idx+1))))
    return z_prime

def locationx(idx):
    #inputs: the index of the element in the matrix
    #returns the z-location of the element on the wing according to the equation given in the assignment
    idx += 1    #add the += 1 as python indexing begins at 0 rather than at 1
    #l_a = 1 #returns as a fraction of chord length. use the lenght of C_a in meters to get a value as such
    theta_x = lambda i : (i-1)*np.pi/41    
    x_prime = 0.5*(l_a/2*(1 - np.cos(theta_x(idx))) + l_a/2*(1 - np.cos(theta_x(idx+1))))
    return x_prime

def MapAeroLoading(filename):
    with open(filename) as aerodynamicloadcrj700:
        #Generate a matrix containing all the data poitns from aerodyanmicloadcrj700
        AeroLoading = np.zeros([81, 41])
        i = 0
        for line in aerodynamicloadcrj700.readlines():
            line = line.strip().split(',')
            AeroLoading[i,:] = np.float_(line)
            i += 1
        return AeroLoading
    
AeroLoading = MapAeroLoading("aerodynamicloadcrj700.dat")
    
def findIndex(loc_x, loc_z):
    #input: x, z; location of a point on the aileron surface
    #outputs: index_x, index_z: the indices of the location in the weight matrix (p11)

    def findz(loc_z):
        for index_z in range(0, 82):
            if locationz(index_z + 1) < loc_z <= locationz(index_z):
                return index_z
        if loc_z > locationz(0):
            return 0
        else:
            return 81
    def findx(loc_x):
        for index_x in range(0, 42):
            if locationx(index_x) <= loc_x < locationx(index_x + 1):
                return index_x
        if loc_x < locationx(0):
            return 0
        else:
            return 41
    return findx(loc_x), findz(loc_z)

def LinearInterpolatePos(Q1, Q2, x_0, x_1, x): 
    #inputs: Q1, Q2, x_0, x_1, x; values of function Q1, Q2 at points x0 and x1, respectively, and a location x where the interpolating function is to be evaluated
    return Q1 + (Q2-Q1)/(x_1 - x_0)*(x-x_0)

def LinearInterpolate(Q1, Q2, x_0, x_1):
    #returns the weights a, b of the interpolating function a + b x
    return np.array([[Q1 - x_0*(Q2 - Q1)/(x_1 - x_0)], [(Q2 - Q1)/(x_1 - x_0)]])



def integrate_1d(x, y, x_max):
    #inputs: x; an array containing all the x-locationx of the points. y; an array containing all the y values of the points. x_i; the value of x to where we integrate.
    #outputs total; the total integrated value until x_i.
    
    #assert len(x) == len(y) #both arrays must have the same size
    
    if x_max > x[-1]:   #if x_max is outside of the range covered by input, we return the total integral of what we can integrate over
        x_max = x[-1]
    if x_max < x[0]:    #if x_max is lower than the lowest x_value, return 0
        return 0
    if len(x) < 2:
        return 0    
    total = 0
    i = 1
    while x[i] < x_max:
        total += (x[i] - x[i-1])*(y[i]+y[i-1])/2
        #print(x[i-1], x[i], total)
        i += 1
    i -= 1
    if x_max > x[i]:
        total += (x_max - x[i])*(LinearInterpolatePos(y[i], y[i+1], x[i], x[i+1], x_max) + y[i])/2
        #print(x[i], x_max, total)
    return total

def integrate_1d_list(x, y, x_max):
    #inputs: x, y; lists containing the locations and values of all data points. x_max; the maximum location until which we integrate
    #outputs: x_new; a list containing all the data locatations up to x_max, and x_max if that is not already in the list. int_list; a list containing the integrated values at each location in x_new
    int_list = []
    x_new = []
    i = 1
    if x_max > x[-1]:
        x_max = x[-1]
    if x_max < x[0]:
        return [0], [x_max]
    if len(x) < 2:
        return [0], [x_max]
    int_list.append(0)
    x_new.append(x[0])
    while x[i] < x_max:
        int_list.append(integrate_1d(x, y, x[i]))
        x_new.append(x[i])
        i += 1
    i -= 1
    if x_max > x[i]:
        int_list.append(integrate_1d(x, y, x_max))
        x_new.append(x_max)

    return int_list, x_new


def integrate_1d_tau(x, y, x_max, x_sc):
    #inputs: x; an array containing all the x-locationx of the points. y; an array containing all the y values of the points. x_i; the value of x to where we integrate.
    #outputs int_x; the total integrated value until x_i.
    
    #assert len(x) == len(y) #both arrays must have the same size
    
    if x_max > x[-1]:   #if x_max is outside of the range covered by input, we return the total integral of what we can integrate over
        x_max = x[-1]
    if x_max < x[0]:    #if x_max is lower than the lowest x_value, return 0
        return 0
    if len(x) < 2:
        return 0    
    total = 0
    i = 1
    while x[i] < x_max:
        total += (x[i] - x[i-1])*(y[i]+y[i-1])/2 * ((x[i] + x[i-1])*0.5 - x_sc)
        print(((x[i] + x[i-1])*0.5 - x_sc))
        #print(x[i-1], x[i], total)
        i += 1
    i -= 1
    if x_max > x[i]:
        total += (x_max - x[i])*(LinearInterpolatePos(y[i], y[i+1], x[i], x[i+1], x_max) + y[i])/2* ((x[i] + x_max)*0.5 - x_sc)
        print((x[i] + x_max)*0.5 - x_sc)
        #print(x[i], x_max, total)
    return total


def integrate_1d_list_tau(x, y, x_max, x_sc):
      #inputs: x, y; lists containing the locations and values of all data points. x_max; the maximum location until which we integrate
    #outputs: x_new; a list containing all the data locatations up to x_max, and x_max if that is not already in the list. int_list; a list containing the integrated values at each location in x_new
    int_list = []
    x_new = []
    i = 1
    if x_max > x[-1]:
        x_max = x[-1]
    if x_max < x[0]:
        return [0], [x_max]
    if len(x) < 2:
        return [0], [x_max]	
    int_list.append(0)
    x_new.append(x[0])
    while x[i] < x_max:
        int_list.append(integrate_1d_tau(x, y, x[i], x_sc))
        x_new.append(x[i])
        i += 1
    i -= 1
    if x_max > x[i]:
        int_list.append(integrate_1d_tau(x, y, x_max, x_sc))
        x_new.append(x_max)

    return int_list, x_new 





# =============================================================================
# Functional Code (that is, not definitions of funcs)
# =============================================================================
    
#Making two arrays with all the x, z locationx
def make_x_and_z():
    x = []
    z = []
    for i in range(0, 41):
        x.append(locationx(i))
    for i in range(0, 81):
        z.append(locationz(i))
    z = z[::-1] #so it goes from -C to 0 instead of 0 to -C
    return x, z

x, z = make_x_and_z()
#Define w_bar
def make_w_bar(AeroLoading = AeroLoading):
    x, z = make_x_and_z()
    w_bar = []
    for index in range(len(x)):
        w_bar = [integrate_1d(z, AeroLoading[:,index], z[-1])] + w_bar
    return w_bar
w_bar = make_w_bar()

def make_tau(x_sc, AeroLoading = AeroLoading):
    x, z = make_x_and_z()
    tau = []
    for index in range(len(x)):
        tau = [integrate_1d_tau(z, AeroLoading[:,index], z[-1], x_sc)] + tau
    return tau



def DoubleIntegral(x_max):
    x = make_x_and_z()[0]
    w_bar = make_w_bar()
    
    return integrate_1d(x, w_bar, x_max)

def ThreeIntegral(x_max):
    x = make_x_and_z()[0]
    # Integration 1
    w_bar = make_w_bar()

    # Integration 2
    int_list_2, x_list_2 = integrate_1d_list(x, w_bar, x_max)

    # Integration 3
    return integrate_1d(int_list_2, x_list_2, x_max)


def FiveIntegral(x_max):
    x = make_x_and_z()[0]
    #Integration 1
    w_bar = make_w_bar()
    
    #Integration 2
    int_list_2, x_list_2 = integrate_1d_list(x, w_bar, x_max)
    
    #Integration 3
    int_list_3, x_list_3 = integrate_1d_list(x_list_2, int_list_2, x_max)
    
    #Integratioan 4
    int_list_4, x_list_4 = integrate_1d_list(x_list_3, int_list_3, x_max)
    
    #Integration 5
    return integrate_1d(x_list_4, int_list_4, x_max)

def DoubleIntegralZSC(x_max, z_sc):
    x = make_x_and_z()[0]
    tau = make_tau(z_sc)
    return integrate_1d(x, tau, x_max)

def TripleIntegralZSC(x_max, z_sc):
    x = make_x_and_z()[0]
    
    #Make Tau (integration 1)
    tau = make_tau(z_sc)
    
    #Integration 2
    int_list_2, x_list_2 = integrate_1d_list(x, tau, x_max)
    
    #Integration 3
    return integrate_1d(x_list_2, int_list_2, x_max)
