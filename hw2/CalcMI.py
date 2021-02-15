import argparse, os, sys
import numpy as np
import pandas as pd
# You can choose to write classes in other python files
# and import them here.


def split_into_bins(data_matrix,bin_num,bin_str):
    gene_numbers,levels = data_matrix.shape
    gene_matrix = []
    if bin_str == 'uniform':
        for i in range(gene_numbers):
            gene = data_matrix[i,:]
            gene_list = []
            max_num = max(gene)
            min_num = min(gene)
            bins = np.linspace(min_num,max_num,bin_num+1)
            for x in gene:
                if x == bins[bin_num]:
                    gene_list.append(bin_num-1)
                elif x == bins[0]:
                    gene_list.append(0)
                else:
                    for j in range(bin_num+1):
                        if x < bins[j]:
                            gene_list.append(j-1)
                            break
            # print('bins',bins)
            # print('gene list', gene_list)
            gene_matrix.append(gene_list)
        return np.array(gene_matrix)
    else:
        gene_matrix = []
        for i in range(gene_numbers):
            gene = data_matrix[i,:]
            index = np.argsort(gene)
            bins_list = np.array_split(index,bin_num)
            gene_list = []
            for j in range(levels):
                for k in range(len(bins_list)):
                    if j in bins_list[k]:
                        gene_list.append(k)
                        break
            gene_matrix.append(gene_list)
            # print(count(gene_list))
        # print(gene_matrix)
        return np.array(gene_matrix)

def construct_count_matrix(bin_num,gene_matrix,gene1_index,gene2_index):
    count_matrix = np.zeros((bin_num,bin_num))
    gene_numbers, levels = gene_matrix.shape
    gene1 = gene_matrix[gene1_index]
    gene2 = gene_matrix[gene2_index]
    for level in range(levels):
        x = gene1[level]
        y = gene2[level]
        count_matrix[x][y] += 1
    count_matrix += 0.1
    return count_matrix

def cal_joint_prob(count_matrix, a,b):
    cur_count = count_matrix[a][b]
    total_count = sum(sum(count_matrix))
    return cur_count/total_count

def cal_single_prob(count_matrix,num,G1=True):
    total_count = sum(sum(count_matrix))
    cur_count = 0
    if G1:
        row = count_matrix[num,:]
        cur_count += sum(row)
    else:
        col = count_matrix[:,num]
        cur_count += sum(col)
    return cur_count/total_count

# This is the main function provided for you.
# Define additional functions to implement mutual information
def main(args):
    # Parse input arguments
    # It is generally good practice to validate the input arguments, e.g.,
    # verify that the input and output filenames are provided
    data_file_path = args.dataset
    bin_num = args.bin_num
    bin_str = args.bin_str
    output_file_path = args.out

    # Where you run your code.
    data = pd.read_table(data_file_path)
    data = np.array(data)
    data = data[:,1:]
    # print(data)
    data_matrix = data.transpose()
    gene_matrix = split_into_bins(data_matrix, bin_num, bin_str)
    print(gene_matrix)

    gene_numbers, levels = gene_matrix.shape

    mutual_inform_list = []
    inform_dict = dict()
    for gene1 in range(gene_numbers):
        for gene2 in range(gene1+1,gene_numbers):
            key_inform = (gene1+1,gene2+1)
            count_matrix = construct_count_matrix(bin_num,gene_matrix,gene1,gene2)
            mutual_inform = 0
            for a in range(bin_num):
                for b in range(bin_num):
                    prob_ab = cal_joint_prob(count_matrix,a,b)
                    prob_a = cal_single_prob(count_matrix,a,True)
                    prob_b = cal_single_prob(count_matrix,b,False)
                    ratio = prob_ab/(prob_a * prob_b)
                    mutual_inform += prob_ab * np.log2(ratio)
            mutual_inform_list.append(mutual_inform)
            inform_dict[key_inform] = mutual_inform
    nparray_mi = np.array(mutual_inform_list)
    indexes = np.argsort(nparray_mi)[::-1]

    keys = list(inform_dict.keys())
    for index in indexes:
        print( str(keys[index]) + ' ' + '%.3f'%(inform_dict[keys[index]]))


# Note: this syntax checks if the Python file is being run as the main program
# and will not execute if the module is imported into a different module
if __name__ == "__main__":
    # Note: this example shows named command line arguments.  See the argparse
    # documentation for positional arguments and other examples.
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('dataset',
                        help='input gene expression data file path',
                        type=str)
    parser.add_argument('--bin_num',
                        help='number of bins',
                        type=int,
                        default=5)
    parser.add_argument('--bin_str',
                        help='binning strategy',
                        type=str,
                        choices={'uniform', 'density'},
                        default='uniform')
    parser.add_argument('--out',
                        help='MI output file path',
                        type=str,
                        default='uniform.txt')

    args = parser.parse_args()
    # Note: this simply calls the main function above, which we could have
    # given any name
    main(args)
