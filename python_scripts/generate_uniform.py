import numpy as np
import argparse


parser = argparse.ArgumentParser()
parser.description = "Generate a sample of a spherical Gaussian distibution in the spatial extent.  Outputs x y z 0 0 0"
parser.formatter_class=argparse.RawDescriptionHelpFormatter
parser.add_argument('R', type=float, 
                    help='The radius of the uniform distribution' +
                    'from which the starting position will be drawn.')
parser.add_argument('N', type=int, 
                    help='The number of macroparticles to draw ' +
                    'from the Gaussian distribution.')
parser.add_argument('--distribution', dest="distribution", type=str, 
                    help='Describes the symmetry of the distribution. Options are spherical or cylindrical.' +
                     'Default is spherical.',default="spherical")
#The tranformed sampling seems to be undersampling the center of the sphere.  I think
#this may be due to round-off error.  So I've re-introduced the longer rejection method.
parser.add_argument('--method', dest="method", type=str, 
                    help='Describes the method used for sampling. Options are rejection or transformed.' +
                     'Default is rejection.',default="rejection")
parser.add_argument('--chirp', dest='chirp', type=float, 
                    help='Applies an exact chirp to the distribution.  Default is 0.',
                    default=0)


args = parser.parse_args()

if args.distribution == "spherical":
  if args.method == "transformed":
    costheta = 2*np.random.uniform(size=args.N) - 1
    theta = np.arccos(costheta)
    phi = 2*np.pi*np.random.uniform(size=args.N)
    u = np.random.uniform(size=args.N)
    r = args.R*np.cbrt(u)
    x = r*np.cos(phi)*np.sin(theta)
    y = r*np.sin(phi)*np.sin(theta)
    z = r*costheta
  elif args.method == "rejection":
    x = np.zeros(args.N)
    y = np.zeros(args.N)
    z = np.zeros(args.N)
    i = 0
    while i < args.N:
      x_samp = np.random.uniform(-1.*args.R,args.R)
      y_samp = np.random.uniform(-1.*args.R,args.R)
      z_samp = np.random.uniform(-1.*args.R,args.R)
      r = np.sqrt(x_samp**2 + y_samp**2 + z_samp**2)
      if r < args.R:
        x[i] = x_samp
        y[i] = y_samp
        z[i] = z_samp
        i += 1
  else:
    raise Exception("That sampling method is unsupoorted")
elif args.distribution == "cylindrical":
  if args.method == "transformed":
    phi = 2*np.pi*np.random.uniform(size=args.N)
    u = np.random.uniform(size=args.N)
    r = args.R*np.sqrt(u)
    x = r*np.cos(phi)
    y = r*np.sin(phi)
    z = np.zeros(x.size)
  if args.method == "rejection":
    x = np.zeros(args.N)
    y = np.zeros(args.N)
    z = np.zeros(args.N)
    i = 0
    while i < args.N:
      x_samp = np.random.uniform(-1.*args.R,args.R)
      y_samp = np.random.uniform(-1.*args.R,args.R)
      r = np.sqrt(x_samp**2 + y_samp**2)
      if r < args.R:
        x[i] = x_samp
        y[i] = y_samp
        i += 1
  else:
    raise Exception("That sampling method is unsupoorted")
else:
  raise Exception("That distribution type is unsupported.")
px = args.chirp*x
py = args.chirp*y
pz = args.chirp*z
for i in range(x.size):
  line = []
  line.append(x[i])
  line.append(y[i])
  line.append(z[i])
  line.append(px[i])
  line.append(py[i])
  line.append(pz[i])
  print " ".join([str(element) for element in line])
