import numpy as np
import sys

def printMatrix(matrix, output_file='out.txt',text='matrix'):
    with open(output_file, "a") as out:
        original_stdout = sys.stdout
        sys.stdout = out
        row,col = matrix.shape
        print(text)
        for i in range(row):
            for j in range(col):
                print("{:.3f}".format(round(matrix[i][j], 3)),end='\t')
            print()
        print()
        sys.stdout = original_stdout

input = np.loadtxt("model2.txt", delimiter='\t')
# print(input)
print(input.shape)

(row,col) = input.shape

printMatrix(input,'part4_out.txt','prob_matrix')

entropy_matrix = np.zeros(input.shape)

for i in range(row):
    for j in range(col):
        entropy_matrix[i,j] = - input[i,j] * np.log2(input[i,j])


printMatrix(entropy_matrix,'part4_out.txt','entropy_matrix')

location_entropy_matrix = np.sum(entropy_matrix,axis=0)
hmax = 1 + np.ones((1,input.shape[1]))


total_height_matrix = hmax - location_entropy_matrix
printMatrix(total_height_matrix,'part4_out.txt','logo total height at each position')

individual_height_matrix = np.zeros(input.shape)

for j in range(col):
    cur_total_height = total_height_matrix[0,j]
    prob_sum = np.sum(input[:,j])
    # print(entropy_sum)
    for i in range(row):
        individual_height_matrix[i,j] = input[i,j]/ prob_sum *cur_total_height

printMatrix(individual_height_matrix,'part4_out.txt','individual character matrix')
