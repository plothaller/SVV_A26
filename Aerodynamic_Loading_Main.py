import numpy as np


# =============================================================================
# Notes:
# file size is 81 x 41
# =============================================================================

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


def interpolateSpline(p11, p12, p21, p22, i_z, i_x):
    #Takes the values at the four corners p11, p21, p12, p22 and returns a function that interpolates the four points.
    #p11 is the location with the lowest x and z value, with p21 having a more negative z value and p12 having a higher x value
    #i_x, i_z are the indices of p11 in the aero loading matrix
    
    x_1 = locationx(i_x)
    x_2 = locationx(i_x+1)
    z_1 = locationz(i_z)
    z_2 = locationz(i_z+1)
    
    A = np.array([[1, x_1, z_1, z_1*x_1],[1, x_1, z_2, x_1*z_2],[1, x_2, z_1, x_2*z_1],[1, x_2, z_2, x_2*z_2]])
    y = np.array([[p11],[p12],[p21],[p22]])
    w = np.linalg.solve(A, y)
    
    return w

        

with open(r"C:\Users\Max van Huffelen\Desktop\Quick Access\University\SVV\aerodynamicloadcrj700.dat") as aerodynamicloadcrj700:
    #Generate a matrix containing all the data poitns from aerodyanmicloadcrj700
    def AeroLoadingMap(file):
        AeroLoading = np.zeros([81, 41])
        i = 0
        for line in aerodynamicloadcrj700.readlines():
            line = line.strip().split(',')
            #print(np.float_(line))
            AeroLoading[i] = np.float_(line)
            i += 1
        return AeroLoading
    AeroLoading = AeroLoadingMap(aerodynamicloadcrj700)
    
    
def generateSplines(AeroLoading = AeroLoading):
    #generates all the splines using interpolateSpline and places them in the final matrix
    WeightMatrix = np.zeros([80, 40, 4])
    for idx_z in range(0, 80):
        for idx_x in range(0, 40):
            p11 = AeroLoading[idx_z, idx_x]
            p12 = AeroLoading[idx_z, idx_x + 1]
            p21 = AeroLoading[idx_z+1, idx_x]
            p22 = AeroLoading[idx_z+1, idx_x+1]
            localSplineWeights = interpolateSpline(p11, p12, p21, p22, idx_z, idx_x)[:].ravel()
            WeightMatrix[idx_z, idx_x, :] = localSplineWeights
    return WeightMatrix

    
    
def findIndex(loc_z, loc_x):
    #input: z, x; location of a point on the aileron surface
    #outputs: index_z, index_x: the indices of the location in the weight matrix (p11)
# =============================================================================
#     def findx(loc_x, guess_index_x):
#         if locationx(guess_index_x+1) >= loc_x >= locationx(guess_index_x):
#             return guess_index_x
#         elif loc_x < locationx(guess_index_x):
#             return findx(loc_x, guess_index_x-1)
#         else:
#             return findx(loc_x, guess_index_x+1)
#     def findz(loc_z, guess_index_z):
#         if locationz(guess_index_z +1) <= loc_z <= locationz(guess_index_z):
#             return guess_index_z
#         elif loc_z > locationz(guess_index_z):
#             return findz(loc_z, guess_index_z -1)
#         else:
#             return findz(loc_z, guess_index_z +1)
#     found_index_z, found_index_x = findz(loc_z, 40), findx(loc_x, 20)
# =============================================================================
    found_index_z = ''
    found_index_x = ''
    for index_z in range(0, 82):
        if locationz(index_z + 1) < loc_z <= locationz(index_z):
            found_index_z = index_z
    for index_x in range(0, 42):
        if locationx(index_x) <= loc_x < locationx(index_x + 1):
            found_index_x = index_x
    if type(found_index_z) == str or type(found_index_x) == str:
        raise AssertionError
    return found_index_z, found_index_x

          
       
WeightMatrix = generateSplines()  
def findInterpolatedValue(loc_z, loc_x, weights = WeightMatrix):
    #inputs: z, x; location of a point on the aileron surface
    #outputs: s; the approximate aerodynamic load at the input location
    
    loc_index_z, loc_index_x = findIndex(loc_z, loc_x)
    localWeights = weights[loc_index_z, loc_index_x, :]
    loc = np.array([1, loc_x, loc_z, loc_x*loc_z])
    
    s = loc@localWeights
    return s
    

def integrate(loc_z_max, loc_x_max, loc_z_min = 0, loc_x_min = 0, weights = WeightMatrix):
    def integrateSpline(index_z, index_x, loc_z_max, loc_x_max, loc_z_min = 0, loc_x_min = 0, weights = WeightMatrix):
        a, b, c, d = weights[index_z, index_x, 0], weights[index_z, index_x, 1], weights[index_z, index_x, 2], weights[index_z, index_x, 3]
        x_min = max(loc_x_min, locationx(index_x))
        x_max = min(loc_x_max, locationx(index_x+1))
        z_min = min(loc_z_min, locationz(index_z))
        z_max = max(loc_z_max, locationz(index_z+1))
        if x_min < x_max and z_min > z_max and locationz(index_z+1) <= z_max <= locationz(index_z) and locationz(index_z+1) <= z_min <= locationz(index_z) and locationx(index_x+1) >= x_max >= locationx(index_x) and locationx(index_x+1) >= x_min >= locationx(index_x):
            Volume = a*(x_max*z_max - x_min*z_min) + b*(x_max*x_max*z_max - x_min*x_min*z_min) + c*(x_max*z_max*z_max - x_min*z_min*z_min) + d*(x_max**2*z_max**2 - x_min**2*z_min**2)
            return Volume
        else:
            raise AssertionError
      
    index_z_min, index_x_min = findIndex(loc_z_min, loc_x_min)
    index_z_max, index_x_max = findIndex(loc_z_max, loc_x_max)
    totalVolume = 0
    for index_z in range(index_z_min, index_z_max+1):
        for index_x in range(index_x_min, index_x_max+1):
            totalVolume += integrateSpline(index_z, index_x, loc_z_max, loc_x_max, loc_z_min, loc_x_min)
    return totalVolume
                

                
                
                
        
 
        