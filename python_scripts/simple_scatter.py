import numpy as np
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(description='Generates a scatter plot of the desired columns from the input file.  Output is the image file.')
parser.add_argument('filename', type=str, help='The filepath to the input file with columns.')
parser.add_argument('-x', dest="x", type=int, help='The column number (beginning with 0) to be used for the x-axis.  Default is 0.', default=0)
parser.add_argument('-y', dest="y", type=int, help='The column number (beginning with 0) to be used for the y-axis.  Default is 1.', default=1)
parser.add_argument('-a', dest="a", type=float, help='The transparency (alpha) of the dots.  Default is 1.', default=1)
parser.add_argument('-d','--delimiter', dest="delimiter", type=str, help='Specifies the delimiter used in the input file.  Default is " ".', default=" ")
parser.add_argument('-o','--output_file', dest="output_file", type=str, help='The path to where the output will be written.  Defualt is simple.png in the current working directory.', default="simple.png")

args = parser.parse_args()

data = np.loadtxt(args.filename,delimiter=args.delimiter).transpose()
plt.scatter(data[args.x],data[args.y],alpha=args.a)
plt.xlim(np.min(data[args.x]),np.max(data[args.x]))
plt.ylim(np.min(data[args.y]),np.max(data[args.y]))

plt.tight_layout()
plt.savefig(args.output_file)
