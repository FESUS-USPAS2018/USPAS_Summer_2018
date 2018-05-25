import numpy as np
import argparse


parser = argparse.ArgumentParser()
parser.description = "Generate a sample of a Gaussian distibution in the spatial extent.  Outputs x y z 0 0 0"
parser.formatter_class=argparse.RawDescriptionHelpFormatter
parser.add_argument('sigma_r', type=float, 
                    help='The radial radius (uniform) or standard deviation (Gaussian) ' +
                    'from which the starting position will be drawn.')
parser.add_argument('sigma_z', type=float, 
                    help='The longitudinal length(uniform) or standard deviation (Gaussian) ' +
                    'from which the starting position will be drawn.')
parser.add_argument('N', type=float, 
                    help='The number of macroparticles to draw ' +
                    'from the Gaussian distribution.')
parser.add_argument('pz', type=float, 
                    help='Longitudinal momentum of the center of the pulse.')
parser.add_argument('--chirp', dest='chirp', type=float, 
                    help='The optional chirp between z and pz (delta_pz/delta_z).', default=None)
parser.add_argument('-z','--zstart', dest='zstart', type=float, 
                    help='The position of the center of the bunch.  Default is 0.', default=0.)
parser.add_argument('-d','--dist', dest='dist', type=str, 
                    help='The distribution from which the position coordinates are drawn.  Default is uniform.  Option is Gaussian.', default="uniform")
                    
args = parser.parse_args()

if args.dist == "uniform":
  phi = 2*np.pi*np.random.uniform(size=int(args.N))
  u = np.random.uniform(size=int(args.N))
  r = args.sigma_r*np.sqrt(u)
  x = r*np.cos(phi)
  y = r*np.sin(phi)
  z = args.sigma_z*np.random.uniform(size=int(args.N))
elif args.dist == "Gaussian":
  cov_matrix = np.identity(3)
  cov_matrix[0][0] = args.sigma_r**2
  cov_matrix[1][1] = args.sigma_r**2
  cov_matrix[2][2] = args.sigma_z**2

  x, y, z = np.random.multivariate_normal([0,0,0], cov_matrix, int(args.N)).transpose()
else:
  raise Exception("Unsuported distribution type.")
px = np.zeros(x.size)
py = np.zeros(y.size)
pz = np.ones(z.size)*args.pz
if args.chirp is not None:
  pz += args.chirp*z
z += args.zstart 
for i in range(x.size):
  line = []
  line.append(x[i])
  line.append(y[i])
  line.append(z[i])
  line.append(px[i])
  line.append(py[i])
  line.append(pz[i])
  print " ".join([str(element) for element in line])
