import numpy as np
import copy
import argparse, os, sys
import math

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
        N += 2
        L = len(line)
    # print(np.array(input_matrix))
    return L,N,np.asarray([np.array(xi) for xi in input_matrix])

def printMatrix(matrix, output_file='out.txt'):
    with open(output_file, "a") as out:
        original_stdout = sys.stdout
        sys.stdout = out

        row,col = matrix.shape
        for i in range(row):
            for j in range(col):
                print(matrix[i][j],end=' ')
            print()
        print()
        sys.stdout = original_stdout

def generate_Prob_Matrix(W,input_matrix,seed):
    prob_matrix = np.zeros((4,W+1))
    np.random.seed(seed)
    row,col = prob_matrix.shape
    for i in range(row):
        for j in range(col):
            prob_matrix[i][j] = np.random.rand()
    return prob_matrix

def get_ACGT_prob(input_matrix):
    prob_A = 0
    prob_C = 0
    prob_G = 0
    prob_T = 0
    count_A = 0
    count_C = 0
    count_G = 0
    count_T = 0
    row,col = input_matrix.shape
    for i in range(row):
        for j in range(col):
            if input_matrix[i,j] == 0:
                count_A += 1
            elif input_matrix[i,j] == 1:
                count_C += 1
            elif input_matrix[i,j] == 2:
                count_G += 1
            else:
                count_T += 1
    total_count = count_A+count_T+count_C+count_G
    prob_A = count_A/total_count
    prob_C = count_C/total_count
    prob_G = count_G/total_count
    prob_T = count_T/total_count
    return total_count,count_A,count_C,count_G,count_T,prob_A,prob_C,prob_G,prob_T

def generate_Z_Matrix(N,W,L,seed):
    np.random.seed(seed)
    z_matrix = np.zeros((N,L-W))
    for i in range(N):
        randNum = np.random.randint(0,L-W,size=1)
        z_matrix[i,randNum] = 1
    return z_matrix

def calZt(input_matrix,prob_matrix,L,W,i):
    Xi = input_matrix[i,:]
    prob_z_list = []
    for start in range(L+1-W):
        probability = 1
        for j in range(L):
            if j in range(0,start):
                letter = Xi[j]
                probability *= prob_matrix[letter][0]
            elif j in range(start,start+W):
                letter = Xi[j]
                probability *= prob_matrix[letter][j-start+1]
            else:
                letter = Xi[j]
                probability *= prob_matrix[letter][0]
        prob_z_list.append(probability)
    sum_prob = sum(prob_z_list)
    norm_prob_z_list = [x / sum_prob for x in prob_z_list]
    return norm_prob_z_list

def get_updated_Z_matrix(input_matrix,prob_matrix,L,W,N):
    updated_z_matrix = []
    for i in range(N):
        norm_prob_z_list = calZt(input_matrix,prob_matrix,L,W,i)
        updated_z_matrix.append(norm_prob_z_list)
    return np.array(updated_z_matrix)

def cal_sum_z_col(z_matrix,input_matrix,L,W,N,letter,pos_in_motif):
    z_score = 0
    for i in range(N):
        Xi = input_matrix[i, :]
        for j in range(L-W):
            # print('pos_in_motif ',pos_in_motif)
            # print('j is ',j)
            if Xi[j + pos_in_motif] == letter:
                # print('find one')
                z_score += z_matrix[i,j]
        # print("sequence is done ----------------", i)
    return z_score

def get_updated_prob_matrix(z_matrix,input_matrix,prob_matrix,L,W,N):
    total_z = np.sum(z_matrix)
    total_count,count_A, count_C,count_G,count_T,prob_A,prob_C,prob_G,prob_T = get_ACGT_prob(input_matrix)
    letter_count_sum = [count_A,count_C,count_G,count_T]
    nck_list = [0,0,0,0]
    nbk = 0
    for letter in range(4):
        count_cur_letter = 0
        cur_letter_sum = letter_count_sum[letter]
        for pos_in_motif in range(W):
            sum_Z_score_at_pos = cal_sum_z_col(z_matrix,input_matrix,L,W,N,letter,pos_in_motif)
            count_cur_letter += sum_Z_score_at_pos
            prob_matrix[letter,pos_in_motif+1] = (sum_Z_score_at_pos + 1)/(total_z + 4)
        nck =cur_letter_sum - count_cur_letter
        nbk += nck
        nck_list[letter] = nck
    for letter in range(4):
        prob_matrix[letter,0] = (nck_list[letter] + 1)  / (nbk + 4)
    return prob_matrix

def detect_threshold(input_matrix, old_z_matrix, new_z_matrix, old_prob_martrix,new_prob_matrix,threshold,W,L,N):
    old_start_list = get_start_list(old_z_matrix)
    new_start_list = get_start_list(new_z_matrix)

    old_log_likelihood = cal_loglikelihood(old_start_list,input_matrix,old_prob_martrix,W,L,N)
    new_log_likelihood = cal_loglikelihood(new_start_list,input_matrix,new_prob_matrix,W,L,N)

    diff = abs(new_log_likelihood-old_log_likelihood)

    if diff < threshold:
        return True
    else:
        return False

def cal_loglikelihood(start_list,input_matrix,prob_matrix,W,L,N):
    log_likelihood = 0
    for i in range(N):
        Xi = input_matrix[i,:]
        start_pos = start_list[i]
        for k in range(start_pos):
            letter = Xi[k]
            cur_prob = prob_matrix[letter,0]
            log_likelihood += math.log2(cur_prob)
        for j in range(W):
            pos = start_pos + j
            letter = Xi[pos]
            cur_prob = prob_matrix[letter,j+1]
            log_likelihood += math.log2(cur_prob)
        for m in range(start_pos + W,L):
            letter = Xi[m]
            cur_prob = prob_matrix[letter,0]
            log_likelihood += math.log2(cur_prob)
    return log_likelihood

def EM(input_matrix,L,W,N,threshold,seed):
    z_matrix= generate_Z_Matrix(N,W,L,seed)
    old_z_matrix = copy.deepcopy(z_matrix)
    prob_matrix = generate_Prob_Matrix(W,input_matrix,seed)
    prob_matrix = get_updated_prob_matrix(z_matrix,input_matrix,prob_matrix,L,W,N)
    old_p = copy.deepcopy(prob_matrix)

    z_matrix = get_updated_Z_matrix(input_matrix, old_p, L, W, N)
    new_z_matrix = copy.deepcopy(z_matrix)

    new_p = get_updated_prob_matrix(z_matrix, input_matrix, prob_matrix, L, W, N)
    count = 0

    while(detect_threshold(input_matrix,old_z_matrix,new_z_matrix,old_p,new_p,threshold,W,L,N) == False):
        old_z_matrix = copy.deepcopy(z_matrix)
        z_matrix = get_updated_Z_matrix(input_matrix, new_p, L, W, N)
        new_z_matrix = copy.deepcopy(z_matrix)

        old_p = copy.deepcopy(new_p)
        new_p = get_updated_prob_matrix(z_matrix, input_matrix, new_p, L, W, N)
        count += 1
    print('em iter: ', count)
    return new_p, z_matrix

# This is the main function provided for you.
# Define additional functions to implement MEME

def findStartPosition(z_matrix):
    start_list = []
    row, col = z_matrix.shape
    for i in range(row):
        seq = z_matrix[i,:]
        order = np.argsort(seq)
        start_list.append(order[-1])
    return start_list

def getMotifSeq(start_list,input_matrix,W,N):
    motif_list = []
    for i in range(N):
        motif = []
        start_position = start_list[i]
        for k in range(W):
            pos = start_position+k
            code = input_matrix[i,pos]
            motif.append(getLetter(code))
        motif_list.append(motif)
    return motif_list

def getLetter(code):
    if code == 0:
        return 'A'
    elif code == 1:
        return 'C'
    elif code == 2:
        return 'G'
    else:
        return 'T'

def printToFile(matrix,output_file):
    with open(output_file, "w") as out:
        original_stdout = sys.stdout
        sys.stdout = out
        row,col = matrix.shape
        for i in range(row):
            print(getLetter(i), end='\t')
            for j in range(col):
                print("{:.3f}".format(round(matrix[i][j], 3)),end='\t')
            print()
        sys.stdout = original_stdout

def printToFile2(motif_list,output_file):
    with open(output_file, "w") as out:
        original_stdout = sys.stdout
        sys.stdout = out
        for pos in motif_list:
            print(pos)
        sys.stdout = original_stdout

def printToFile3(motif_list,input_matrix, N, L, W, output_file):
    with open(output_file, "w") as out:
        original_stdout = sys.stdout
        sys.stdout = out
        motifs = []
        for i in range(N):
            start_pos = motif_list[i]
            seq = input_matrix[i]
            motif = ''
            for j in range(W):
                cur_pos = start_pos + j
                motif += getLetter(seq[cur_pos])
            motifs.append(motif)
        for motif in motifs:
            print(motif)
        sys.stdout = original_stdout

def get_start_list(z_matrix):
    N,col = z_matrix.shape
    start_list = []
    for i in range(N):
        cur_seq = z_matrix[i]
        sorted_seq = np.argsort(cur_seq)
        start_list.append(sorted_seq[-1])
    return start_list

def get_subsequence (input_matrix,W,N,L):
    subsequences = []
    for i in range(N):
        cur_seq = input_matrix[i]
        for j in range(0,L-W+1):
            subsequence = []
            for k in range(0,W):
                letter = cur_seq[j+k]
                subsequence.append(letter)
            subsequences.append(subsequence)
    return subsequences

def main(args):
    seq_file_path = args.sequences_filename
    W = args.width
    model_file_path = args.model
    position_file_path = args.positions
    subseq_file_path = args.subseqs
    L, N, input_matrix = getSeq(seq_file_path)
    seed_list = []
    likelihood_list = []
    for iter in range(100):
        cur_seed = iter * 23 + 13
        print('loop ', iter)
        seed_list.append(cur_seed)
        final_p_matrix, final_z_matrix = EM(input_matrix,L,W,N,1e-3,cur_seed)
        start_list = findStartPosition(final_z_matrix)
        likelihood = cal_loglikelihood(start_list, input_matrix, final_p_matrix, W, L, N)
        likelihood_list.append(likelihood)
    likelihood_arr = np.array(likelihood_list)
    sort_index = np.argsort(likelihood_arr)
    best_seed = seed_list[sort_index[-1]]
    best_likelihood = likelihood_arr[sort_index[-1]]
    # print(best_seed)
    print(best_likelihood)
    np.random.seed(best_seed)
    final_p_matrix, final_z_matrix = EM(input_matrix, L, W, N, 1e-8, best_seed)
    start_list = findStartPosition(final_z_matrix)

    printToFile(final_p_matrix,model_file_path)
    printToFile2 (start_list,position_file_path)
    printToFile3 (start_list, input_matrix, N, L, W, subseq_file_path)

    best_likelihood = cal_loglikelihood(start_list, input_matrix, final_p_matrix, W, L, N)
    print('--------------best likelihood---------------------------------')
    print(best_likelihood)


# Note: this syntax checks if the Python file is being run as the main program
# and will not execute if the module is imported into a different module
if __name__ == "__main__":
    # Note: this example shows named command line arguments.  See the argparse
    # documentation for positional arguments and other examples.
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('sequences_filename',
                        help='sequences file path.',
                        type=str)
    parser.add_argument('--width',
                        help='width of the motif.',
                        type=int,
                        default=6)
    parser.add_argument('--model',
                        help='model output file path.',
                        type=str,
                        default='model.txt')
    parser.add_argument('--positions',
                        help='position output file path.',
                        type=str,
                        default='positions.txt')
    parser.add_argument('--subseqs',
                        help='subsequence output file path.',
                        type=str,
                        default='subseqs.txt')

    args = parser.parse_args()
    # Note: this simply calls the main function above, which we could have
    # given any name
    main(args)
