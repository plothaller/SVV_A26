import numpy as np


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

with open(r"C:\Users\Max van Huffelen\Desktop\Quick Access\University\SVV\aerodynamicloadcrj700.dat") as aerodynamicloadcrj700:
    #Generate a matrix containing all the data poitns from aerodyanmicloadcrj700
    def MapAeroLoading(file):
        AeroLoading = np.zeros([81, 41])
        i = 0
        for line in aerodynamicloadcrj700.readlines():
            line = line.strip().split(',')
            AeroLoading[i] = np.float_(line)
            i += 1
        return AeroLoading
    AeroLoading = MapAeroLoading(aerodynamicloadcrj700)
    
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


def BilinearInterpolate(Q11, Q21, Q12, Q22, x_0, z_0, x_1, z_1, x, z, z_negative = True):
    #returns the approximate value at a location (x,z) between points Qab at (x_a, y_b)
    # if (x,z) lies outside [(x_0, z_0), (x_1, z_1)] the nearest point inside the domain will be taken instead
    x = max(min(x_1, x), x_0)
    if z_negative:
        z = min(max(z_1, z), z_0)
    else:
        z = max(min(z_1, z), z_0)
    return (Q11*(x_1 - x)*(z_1 - z) + Q21*(x - x_0)*(z_1 - z) + Q12*(x_1 - x)*(z - z_0) + Q22*(x - x_0)*(z - z_0))/((x_1 - x_0)*(z_1 - z_0))
    
    

def integrate(x_0, z_0, x_1, z_1, AeroLoad = AeroLoading):
    #Integrate from (x_0, z_0) to (x_1, z_1)
    index_x_0, index_z_0 = findIndex(x_0, z_0)
    index_x_1, index_z_1 = findIndex(x_1, z_1)
    
    def integrateComponent(index_x, index_z, x_0 = x_0, z_0 = z_0, x_1 = x_1, z_1 = z_1, AeroLoad = AeroLoad):
        #define p11, p12, p21, p22 where p11 has lowest x, z indices, p12 has higher z index, p21 has higher x index, p22 has higher x and z index
        
        #Define the loads at the cornerpoints of the component section
        Q11 = AeroLoad[index_z, index_x]
        Q12 = AeroLoad[index_z+1, index_x]
        Q21 = AeroLoad[index_z, index_x+1]
        Q22 = AeroLoad[index_z+1,index_x+1]
        #Define the locations of the cornerpoints
        x_lower = locationx(index_x)
        x_upper = locationx(index_x+1)
        z_lower = locationz(index_z)    #this is lower INDEX, not lower value
        z_upper = locationz(index_z+1)
        #Find the values of aero loading (estimated) at each of the corners of the SUB-AREA WE ARE INTERPOLATING, (not the section of aero data)
        p11 = BilinearInterpolate(Q11, Q21, Q12, Q22, x_lower, z_lower, x_upper, z_upper, x_0, z_0)
        p12 = BilinearInterpolate(Q11, Q21, Q12, Q22, x_lower, z_lower, x_upper, z_upper, x_0, z_1)
        p21 = BilinearInterpolate(Q11, Q21, Q12, Q22, x_lower, z_lower, x_upper, z_upper, x_1, z_0)
        p22 = BilinearInterpolate(Q11, Q21, Q12, Q22, x_lower, z_lower, x_upper, z_upper, x_1, z_1)
        #integrate twice using the trapezoid rule
        return (x_upper-x_lower)*(z_upper-z_lower)*(p11 + p12 + p21 + p22)/4
    
    totalIntegrated = 0
    
    for index_x in range(min(index_x_0, index_x_1), min(max(index_x_1, index_x_0), 40)):
        for index_z in range(min(index_z_0, index_z_1), min(max(index_z_0, index_z_1), 80)):
            #add the integrated value of one component to the total. Do this for each component to get the total integrated value
            totalIntegrated += integrateComponent(index_x, index_z)
    return totalIntegrated * np.sign((index_x_1-index_x_0)*(index_z_1 - index_z_0))

















