import numpy as np
import sys

def getSeq(filename):
    file = open(filename, 'r')
    Lines = file.readlines()
    input_matrix = []
    N = 0
    L = 0
    for line in Lines:
        line = line[:-1]
        num_list = processSeq(line)
        input_matrix.append(num_list)
        N += 1
        L = len(line)
    return L,N,np.array(input_matrix)

def processSeq(x):
    x_upper= x.upper()
    num_list = []
    for i in x_upper:
        if i == 'A':
            num_list.append(0)
        elif i == 'C':
            num_list.append(1)
        elif i == 'G':
            num_list.append(2)
        else:
            num_list.append(3)
    return num_list

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

def generate_Prob_Matrix(W):
    prob_matrix = np.zeros((4,W+1))
    row,col = prob_matrix.shape
    for i in range(row):
        for j in range(col):
            prob_matrix[i][j] = 0.25
    return prob_matrix

def generate_n_matrix(input_matrix,W,N,L,exclude_seq,start_list):
    n_matrix = np.zeros((4, W + 1))
    for i in range(N):
        if i == exclude_seq:
            continue

        cur_start = start_list[i]
        for j in range(cur_start):
            letter = input_matrix[i][j]
            n_matrix[letter][0] += 1
        for j in range(W):
            cur_pos = cur_start + j
            letter = input_matrix[i][cur_pos]
            n_matrix[letter][j+1] += 1
        for j in range(cur_start+W,L):
            letter = input_matrix[i][j]
            n_matrix[letter][0] += 1
    return n_matrix

def update_prob_matrix(n_matrix,prob_matrix,N,W,L):
    for letter in range(4):
        prob_matrix[letter][0] = (n_matrix[letter][0] + 1) / ((N-1)*(L-W) + 4 )
        for j in range(1,W+1):
            prob_matrix[letter][j] = (n_matrix[letter][j] + 1 )/ (N-1+4)
    return prob_matrix

def cal_LR(input_matrix,start,W,prob_matrix,i):
    prob_nom = 1
    prob_dom = 1
    for j in range(W):
        prob_pos = j + 1
        seq_pos = start + j
        letter = input_matrix[i][seq_pos]
        prob_nom *= prob_matrix[letter][prob_pos]
        prob_dom *= prob_matrix[letter][0]
    return prob_nom/prob_dom

def get_LR_sum(input_matrix,prob_matrix,L,W,N,i):
    total_LR = 0
    for k in range(L-W+1):
        total_LR += cal_LR(input_matrix,k,W,prob_matrix,i)
    return total_LR

def update_prob_palindrome(n_matrix,prob_matrix,N,W,L):
    dom_num = 0
    for letter in range(4):
        dom_num += n_matrix[letter][1]
    dom_num = (dom_num + 4)*2

    for letter in range(4):
        op_letter = 3-letter
        for j in range(1,int(W/2)):
            op_j = W + 1 - j
            nom_num = n_matrix[letter][j] + 1 + n_matrix[op_letter][op_j] + 1
            prob_matrix[letter][j] = nom_num/dom_num
            prob_matrix[op_letter][op_j] = nom_num/dom_num
    return prob_matrix

L,N,input_matrix = getSeq('seq.txt')
W = 4
exclude_seq = 4 #seq start from 0.... N-1
start_list = [1,3,1,0,2,1,0,1,2,3]
prob_matrix = generate_Prob_Matrix(4)
printMatrix(input_matrix,'part2_out.txt','sequence matrix')
printMatrix(prob_matrix,'part2_out.txt','initialized prob matrix')
n_matrix = generate_n_matrix(input_matrix,W,N,L,exclude_seq,start_list)
printMatrix(n_matrix,'part2_out.txt','N matrix')
prob_matrix = update_prob_matrix(n_matrix,prob_matrix,N,W,L)
printMatrix(prob_matrix,'part2_out.txt','updated probability matrix')

LR = cal_LR(input_matrix,1,4,prob_matrix,4)
total_LR = get_LR_sum(input_matrix,prob_matrix,L,W,N,4)
# print(total_LR)
prob_start = LR/total_LR
with open('part2_out.txt', "a") as out:
    original_stdout = sys.stdout
    sys.stdout = out
    print('probability of choosing ai=2:')
    print("{:.3f}".format(round(prob_start, 3)))
    print()
    sys.stdout = original_stdout

palidrome_prob_matrix = update_prob_palindrome(n_matrix,prob_matrix,N,W,L)
printMatrix(palidrome_prob_matrix,'part2_out.txt','updated palidrome probability matrix')
