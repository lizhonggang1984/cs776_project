# Note: this is the Python syntax for importing modules
# We sometimes use package as a synonym for module, and formally
# a package is a module with a path on the file system.
# Modules may have submodules, which are delimited with .
# i.e. module.submodule
import argparse, os, sys

# Note: 'as' provides an alias for a module so that you can
# access its functions with syntax like np.mean() instead of numpy.mean()
import numpy as np

np.random.seed(776)

# Note: this is the syntax for Python functions. You do not need to
# specify a return type in advance. Multiple values can be returned
# using tuples such as (value1, value2, value3).
def read_file(input_filename='data.txt'):
    # Note: informally, 'with' statements are used for file reading and
    # writing, among other things. They guarantee that the file is properly
    # closed regardless of whether the code block runs normally or throws
    # an exception. The line below opens the file in read mode and creates
    # a file object.
    with open(input_filename, "r") as in_file:
        # Read the file contents.
        # For normal NumPy data, you can load by np.load.
        # But for practice, you need to read them from plain text here.
        data = np.loadtxt(in_file)
    return data

def hw0_test(input_file, output_file):
    # Note: out is the output file you can write to with the print function.
    with open(output_file, "w") as out:
        original_stdout = sys.stdout
        sys.stdout = out
        #####
        #Part 1: strings
        #####
        # Create variables first_name and last_name with your name.
        first_name = 'Zhonggang'
        last_name = 'Li'
        # Create full_name by concatenating the variables with a tab
        full_name = first_name + "  " + last_name

        # separating them and print it.
        full_name_list = full_name.split()

        print('First name: ',full_name_list[0])
        print()
        print('Last name: ',full_name_list[1])
        print()
        # Transform the full_name string to all uppercase letters and
        # print it.
        upper_full_name = full_name.upper()
        print('Full name: ', upper_full_name)
        print()

        #####
        # Part 2: lists
        #####
        # Initialize a list x with four zeros.
        x = [0 for i in range(4)]
        # Prepend one 1 to the head and append one 1 to the tail. Print x.
        x = [1] + x
        x.append(1)
        print('List of x: ', x)
        print()
        # Set y to be the first four elements of x.
        y = x[:4]
        # Change the second to last element of y to 2. Print y.
        for i in range(len(y)):
            if i >= 1:
                y[i] = 2
        print('List of y:  ', y)
        print()
        # Calculate the product of the elements in y.
        # Pass (skip over) the element if it is 0.  Print the result.

        product = 1
        for i in range(len(y)):
            if y[i] != 0:
                product *= y[i]
        print('Product is %d' %product)
        print()

        # Note: Python strings can be indexed in the same manner as lists.
        course_str = "Advanced Bioinformatics"
        # Print the index of 'B'.
        index_of_B = course_str.index('B')
        print('Index for B is %d' %index_of_B)
        print()
        # Print the number of occurrences of "i" and assign it to variable k.
        k = 0
        for i in range(len(course_str)):
            if course_str[i] == 'i':
                k += 1
        print('i occured %d times in course_str'%k)
        print()
        # Print the substring of length k starting with 'B'.
        sub_str = course_str[index_of_B:index_of_B+k]
        print('String starts with B: ', sub_str)
        print()
        #####
        # Part 3: dictionaries and sets
        #####
        # Note: This is set syntax.
        keys = {"a", "b", "c", "d"}
        # Create a dictionary called hash_map.
        # Map the chararacters a-d to 1-4. Save the mapping in hash_map.
        hash_map = dict()
        hash_map['a'] = 1
        hash_map['b'] = 2
        hash_map['c'] = 3
        hash_map['d'] = 4

        # Check if "e" exists in hash_map. If not, map it to 5.
        if 'e' in hash_map:
            print('e exists in hash map')
            print()
        else:
            hash_map['e'] = 5
        # Print all key-value pairs in format <key:value> like "a:1".
        print(hash_map)
        print()
        # Map "e" to 6.
        hash_map['e'] = 6
        # Print all key-s
        print(hash_map)
        print()
        #####
        # Part 4: NumPy arrays
        # 
        # Note: you may write a function print_matrix(matrix, output_file) 
        # to print matrices in the desired format.
        #####
        u = [1, 2, -3]
        v = [3, -2, 1]
        # Convert u and v into NumPy arrays.

        u = np.array(u)
        v = np.array(v)

        # Calculate a, the inner product of u and v. Print a.
        # Hint: a is a scalar.
        a = u.dot(v)
        print('innter product a is %d'%a)
        print()

        # Calculate B, the outer product of u and v. Print B.
        # Hint: B is a 3 by 3 matrix.
        B = np.outer(u,v)
        print('outer product B is: ')
        printMatrix(B)
        print()

        # Calculate C = a * B. Print C.
        print('scaler with matrix product C is ')
        C = a*B
        printMatrix(C)
        print()

        # Create R, a 3 by 3 random matrix. Print R.
        R = np.random.rand(3,3)
        print('R is')
        printMatrix(R)
        print()

        # Calculate the matrix product RC. Print the result.
        RC = R.dot(C)
        print('RC is')
        printMatrix(RC)
        print()

        # Calculate the matrix product CR. Print the result.
        CR = C.dot(R)
        print('CR is')
        printMatrix(CR)
        print()

        # Calculate the elementwise product of R and C. Print the result.
        ElementWise = np.multiply(R,C)
        print('Element product is ')
        printMatrix(ElementWise)
        print()

        #####
        # Part 5: NumPy sorting and binning
        #####
        # Complete the read_file function and call it to
        # read the matrix D from input_file into a NumPy array.
        data = read_file()

        # Sort the first column of D in ascending order.
        # Print the first three values of the sorted list.

        first_col = data[:,0]
        argindexes = first_col.argsort()
        sorted_data = []
        for argindex in argindexes:
            sorted_data.append(first_col[argindex])
        sorted_data = np.array(sorted_data)
        print('first_three elements are: ')
        for i in range(3):
            print(sorted_data[i], end=' ')
        print()

        # import random as rd
        # data = [rd.randint(a=100, b=1000) for _ in range(20)]
        # bins = [200, 300, 400, 500, 600, 700, 800, 900, 1000]
        # print('data:', data)
        # print('bins:', bins)
        # print('np.digitize(data,bins):', np.digitize(data, bins))

        # Assign the sorted values into five bins with equal width.
        # The left edge of the first bin is the min value.
        # The right edge of the last bin is the max value.
        # For each bin, print
        #  1) its left and right edges, and
        #  2) the values that fall inside the bin.

        max_num = max(sorted_data)
        min_num = min(sorted_data)

        bins = np.linspace(min_num,max_num,6)
        print()
        x = sorted_data
        inds = np.digitize(x,bins)

        print('bin sorting with equal bin width: ')
        left_set = set()
        right_set = set()
        for ccc in range(x.size-1):
            left_set.add(bins[inds[ccc]-1])
            right_set.add(bins[inds[ccc]])
        left_list = list(left_set)
        right_list = list(right_set)
        left_list.sort()
        right_list.sort()
        left_list.append(bins[5])
        right_list.append(bins[5])
        bin_list = [[],[],[],[],[]]

        for n in range(x.size):
            cur_num = x[n]
            for i in range(5):
                if cur_num <= right_list[i]:
                    bin_list[i].append(cur_num)
                    break

        for i in range(5):
            print('bin %d values: '%i, end=' ')
            print(bin_list[i])
            print('left edge: %.4f; right edge: %.4f' %(left_list[i],right_list[i]))
        # Assign the sorted values into five bins such that each bin
        # contains the same number of elements.
        # Print the values that fall inside each bin.
        print()
        a,b,c,d,e = np.split(sorted_data,5)
        print('bin sorting with equal numbers in each bin: ')
        print(a)
        print('left edge: %.4f; right edge: %.4f' % (a[0], a[-1]))
        print(b)
        print('left edge: %.4f; right edge: %.4f' % (b[0], b[-1]))
        print(c)
        print('left edge: %.4f; right edge: %.4f' % (c[0], c[-1]))
        print(d)
        print('left edge: %.4f; right edge: %.4f' % (d[0], d[-1]))
        print(e)
        print('left edge: %.4f; right edge: %.4f' % (e[0], e[-1]))
        print()

        # Sort the last row of D in descending order.
        # Print the last three values of the sorted list.
        rev_sorted_data = sorted_data[::-1]
        print('first_three elements descending are: ')
        for i in range(3):
            print(rev_sorted_data[i], end=' ')
        print()

        sys.stdout = original_stdout

def printMatrix(matrix, output_file='out.txt'):
    with open(output_file, "w") as out:
        row,col = matrix.shape
        for i in range(row):
            for j in range(col):
                print(matrix[i][j],end=' ')
            print()

def main(args):
    input_file = args.inputfile
    output_file = args.outputfile
    hw0_test(input_file, output_file)


# Note: this syntax checks if the Python file is run as the main program
# and will not execute if the module is imported into a different module.
if __name__ == "__main__":
    # Note: this example shows named command line arguments. See the argparse
    # documentation for positional arguments and other examples.
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    # Note: you can use ' or " for strings.
    parser.add_argument('--inputfile',
                        help='input data file path.',
                        type=str,
                        default='data.txt')
    parser.add_argument('--outputfile',
                        help='output file path.',
                        type=str,
                        default='out.txt')

    args = parser.parse_args()
    # Note: this simply calls the main function above, which we could have
    # given any name.
    main(args)
