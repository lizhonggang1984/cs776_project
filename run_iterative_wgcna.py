#!/usr/bin/env python
import os
'''Convenience wrapper for running iterativeWGCNA directly from source tree.'''
# import readline

from iterativeWGCNA.cmlargs import parse_command_line_args
from iterativeWGCNA.iterativeWGCNA import IterativeWGCNA

# if __name__ == '__main__':
cmlArgs = parse_command_line_args()
alg = IterativeWGCNA(cmlArgs)
alg.run()



