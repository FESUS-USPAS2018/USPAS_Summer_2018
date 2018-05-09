import numpy as np
from scipy import constants
import argparse

def calc_cov_matrix(x,p):
  xmean = np.mean(x)
  pmean = np.mean(p)
  xsqmean = np.mean(x*x)
  psqmean = np.mean(p*p)
  xpmean = np.mean(x*p)
  
  varx = xsqmean - xmean**2
  varp = psqmean - pmean**2
  covxp = xpmean - xmean*pmean

  return np.array([[varx,covxp],[covxp,varp]])

def calc_emittance(cov_mat):
  return (cov_mat[0][0]*cov_mat[1][1] - cov_mat[0][1]**2)

parser = argparse.ArgumentParser(description='Prints the emittance calculated from the input file.  Output is a number to the terminal.')
parser.add_argument('filepath', type=str, help='The filepath to phase volume.  The input file should have ","-deliminated columns x,y,z,px,py,pz.  Units are assumed to be mks.')
#parser.add_argument('NperM', type=float, help='Number of electrons per macroparticle.')



m = constants.physical_constants["electron mass"][0]
c = constants.physical_constants['speed of light in vacuum'][0]

args = parser.parse_args()


z, p = np.loadtxt(open(args.filepath,'r'),delimiter=',').transpose()
p = p/(m*c)
cov_mat = calc_cov_matrix(z,p)
emittance = np.sqrt(calc_emittance(cov_mat))
print str(emittance)
