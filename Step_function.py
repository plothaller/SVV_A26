import numpy as np


def step_matrix(x, step_list):
	matrix = np.identity(len(step_list))
	for i in range(0,len(step_list)):
		if step_list[i] > x:
			matrix[i] = np.zeros(shape=(len(step_list)))
	return matrix

print(step_matrix(5,(1,2,8,5)))

	
	
