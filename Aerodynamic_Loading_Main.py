import numpy as np
import math


# =============================================================================
# Notes:
# file size is 81 x 41
#
# =============================================================================

with open(r"C:\Users\Max van Huffelen\Desktop\Quick Access\University\SVV\aerodynamicloadcrj700.dat") as aerodynamicloadcrj700:
    i = 0
    matrixAeroLoading = np.zeros([81, 41])
    for line in aerodynamicloadcrj700.readlines():
        matrixAeroLoading[i,:] = line.strip().split(',')
        i += 1
        if i == 1:
            print(line.strip().split(','))
    print(matrixAeroLoading[0,:])
 
        