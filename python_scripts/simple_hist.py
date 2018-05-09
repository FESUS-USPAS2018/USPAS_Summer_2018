import numpy as np
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(description='Generates a histogram plot of the desired column from the input file.  The output is the image file.')
parser.add_argument('filename', type=str, help='The filepath to the input file with columns.')
parser.add_argument('-x', dest="x", type=int, help='The column number (beginning with 0) to be used for the x-axis.  Default is 0.', default=0)
parser.add_argument('-b', dest="bins", type=int, help='The number of bins.  Default is 20.', default=20)
parser.add_argument('-d','--delimiter', dest="delimiter", type=str, help='Specifies the delimiter used in the input file.  Default is " ".', default=" ")
parser.add_argument('-o','--output_file', dest="output_file", type=str, help='The path to where the output will be written.  Defualt is simple.png in the current working directory.', default="simple.png")

args = parser.parse_args()

data = np.loadtxt(args.filename,delimiter=args.delimiter).transpose()
plt.hist(data[args.x],bins=args.bins)

plt.tight_layout()
plt.savefig(args.output_file)
