import numpy as np
import sys
from scipy.special import gamma, factorial

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
    prob_matrix = np.zeros((4, W + 1))
    row, col = prob_matrix.shape
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

def import_alpha_matrix(filename):
    file = open(filename, 'r')
    Lines = file.readlines()
    alpha_matrix = []
    for line in Lines:
        line = line[:-1]
        row_list = []
        for alpha in line:
            if alpha ==' ':
                continue
            row_list.append(int(alpha))
        alpha_matrix.append(row_list)
    return np.array(alpha_matrix)

def cal_p_nk_in_alpha(n_matrix,alpha_matrix,alpha_upper_index,k):
    p_nk_with_alpha_upper = 1
    n_k = sum (n_matrix[:,k])
    alpha_sum = sum(alpha_matrix[:,alpha_upper_index])
    # print(n_k)
    p_nk_with_alpha_upper *= (gamma(n_k + 1) * gamma(alpha_sum) / gamma(n_k + alpha_sum))

    for c in range(4):
        n_c_k = n_matrix[c,k]
        alpha_c = alpha_matrix[c,alpha_upper_index]
        p_nk_with_alpha_upper *= (gamma(n_c_k + alpha_c) / (gamma(n_c_k+1) * gamma(alpha_c)))
    return p_nk_with_alpha_upper


def get_normalized_p_alpha_in_nk(n_matrix,alpha_matrix,k):
    _, alpha_num = alpha_matrix.shape
    p_alpha_in_nk_list = []
    total_p = 0
    for j in range(alpha_num):
        total_p += cal_p_nk_in_alpha(n_matrix,alpha_matrix,j,k)
    for j in range(alpha_num):
        p_alpha_j_in_nk = cal_p_nk_in_alpha(n_matrix,alpha_matrix,j,k)/total_p
        p_alpha_in_nk_list.append(p_alpha_j_in_nk)
    # print('p_alpha_list ',p_alpha_in_nk_list)
    return p_alpha_in_nk_list

def generate_d_matrix(n_matrix,alpha_matrix,W):
    d_matrix = np.zeros((4,W+1))
    _, alpha_num = alpha_matrix.shape
    for c in range(4):
        for k in range(W+1):
            d = 0
            p_alpha_in_nk_list = get_normalized_p_alpha_in_nk(n_matrix,alpha_matrix,k)
            for upper in range(alpha_num):
                d += (p_alpha_in_nk_list[upper] * alpha_matrix[c,upper])
            d_matrix[c,k] = d
    #print(d_matrix)
    return d_matrix

def update_prob_matrix_with_d_matrix(n_matrix,prob_matrix,d_matrix,N,W,L):
    for k in range(W+1):
        total_nbk = 0
        for c in range(4):
            total_nbk += (n_matrix[c,k] + d_matrix[c,k])
        for c in range(4):
            p_ck = (n_matrix[c,k] + d_matrix[c,k]) / total_nbk
            prob_matrix[c,k] = p_ck
    return prob_matrix


L,N,input_matrix = getSeq('seq.txt')
W = 4
exclude_seq = 4 #seq start from 0.... N-1
start_list = [1,3,1,0,2,1,0,1,2,3]
prob_matrix = generate_Prob_Matrix(4)
# printMatrix(input_matrix,'part3_out.txt','sequence matrix')
# printMatrix(prob_matrix,'part3_out.txt','initialized prob matrix')
n_matrix = generate_n_matrix(input_matrix,W,N,L,exclude_seq,start_list)
printMatrix(n_matrix,'part3_out.txt','N matrix')

#print(n_matrix)
alpha_matrix = import_alpha_matrix('alpha_matrix.txt')
#print(alpha_matrix)

# p_nk_with_alpha_upper_1 = cal_p_nk_in_alpha(n_matrix,alpha_matrix,0,0)
# print(p_nk_with_alpha_upper_1)
d_matrix = generate_d_matrix(n_matrix,alpha_matrix,W)
printMatrix(d_matrix,'part3_out.txt','d matrix')

updated_prob_matrix = update_prob_matrix_with_d_matrix(n_matrix,prob_matrix,d_matrix,N,W,L)

printMatrix(updated_prob_matrix,'part3_out.txt','updated probablilty PWM matrix')

