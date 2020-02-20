import numpy as np


# =============================================================================
# Notes:
# file size is 81 x 41
# =============================================================================
#global variables C_a and l_a for the chord and length of the aileron, respectively
C_a = 1
l_a = 1



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
    WeightMatrix = generateSplines()
    
    
    def findIndex(loc_z, loc_x):
        #input: z, x; location of a point on the aileron surface
        #outputs: index_z, index_x: the indices of the location in the weight matrix (p11)
        guess_index_z = np.arccos(1 + 2*loc_z/C_a)*81/np.pi + 1 - 0.5
        guess_index_x = np.arccos(1 - 2*loc_x/l_a)*41/np.pi + 1 - 0.5

        
        def validate_guess_z(guess_index_z, z):
            print('guess: ', guess_index_z, '\n guess location: ', locationz(int(guess_index_z)), '\n desired location: ', loc_z)
            guess_index_z = int(guess_index_z)
            if locationz(guess_index_z) <= z <= locationz(guess_index_z +1):
                return guess_index_z
            elif locationz(guess_index_z) > z:
                validate_guess_z(guess_index_z-1,z)
            else:
                validate_guess_z(guess_index_z+1, z)
                
        def validate_guess_x(guess_index_x, x):
            print('guess: ', guess_index_x, '\n guess location: ', locationz(int(guess_index_x)), '\n desired location: ', loc_x)
            guess_index_x = int(guess_index_x)
            if locationz(guess_index_x) <= x <= locationz(guess_index_x +1):
                return guess_index_x
            elif locationz(guess_index_x) > x:
                validate_guess_x(guess_index_x-1,x)
            else:
                validate_guess_x(guess_index_x+1, x)
        index_z, index_x = validate_guess_z(guess_index_z, loc_z), validate_guess_x(guess_index_x, loc_x)
        return index_z, index_x
                
    def findInterpolatedValue(loc_z, loc_x, weights = WeightMatrix):
        #inputs: z, x; location of a point on the aileron surface
        #outputs: s; the approximate aerodynamic load at the input location
        
        #[find at which spline the location is]
        localWeights = weights[loc_index_z, loc_index_x, :]
        loc = np.array([1, loc_x, loc_z, loc_x*loc_z])
        
        s = loc@localWeights
        return s
        
                
    
                
                
                
                
        
 
        