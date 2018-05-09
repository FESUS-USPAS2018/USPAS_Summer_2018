import argparse
import numpy as np
from scipy import constants

def uniform_disc_sample(n,R):
  """
  Returns the x and y positions of points distributed uniformally in 2D.
  args:
    n: the number of positions desired.
    R: the raidus of the desired distribution.  For uniform, this is the radius, for Gaussian, the standard deviation.
  returns:
    x,y : Two np.arrays() that are distributed uniformally in a circle on the x-y plane.
  """
  a = np.random.uniform(size=n)
  b = np.random.uniform(size=n)
  c = np.minimum(a,b)
  d = np.maximum(a,b)
  return (d*R*np.cos(2*np.pi*c/d), d*R*np.sin(2*np.pi*c/d))

def spatial_transverse_distribution(n,R,d="uniform",**kwargs):
  """
  Generates random x,y positions according to the specified distribution.
  args:
    n: the number of positions desired.
    R: the raidus of the desired distribution.  For uniform, this is the radius, for Gaussian, the standard deviation.
    d: the distribution.  Supports 'uniform' and 'Gaussian'.  Default is 'uniform'.
  returns:
    x,y : Two np.arrays that are distributed appropriately in a circle on the x-y plane.
  """
  if d == "uniform":
    x, y = uniform_disc_sample(n,R)
  elif d == "Gaussian":
    x, y = np.random.multivariate_normal([0,0], R**2*np.identity(2), n).transpose()
  else:
    raise Exception("Distribution not supported.")
  return x, y


def photoemission_p(n,E_fermi=5.000E-6,E_work=4.45E-6,E_photon=4.65E-6,e_per_mpart=100,units="mks",**kwargs):
  """
  Simulates the photoemission process n times resulting in n points in momentum spaces assuming non-relativistic dynamics.
  args:
    n: the number of momenta desired.
    E_fermi: the fermi energy measured in MeV.
    E_work: the work function measured in MeV.
    E_photon: the energy of the incoming photon measured in MeV.
    e_per_mpart: the number of electrons per macroparticle.
    units: specifies the desired units for the momenta.  Default is mks (kg*m/s).  Supports mks and MeV (MeV/c).
  returns:
    px, py, pz:  Momentum vetors of length n with the simulated momenta in mks units.
  """
  m_e = constants.physical_constants["electron mass energy equivalent in MeV"][0]
  c = constants.physical_constants["speed of light in vacuum"][0]
  if units == "mks":
    m_mpart = e_per_mpart*constants.physical_constants["electron mass"][0]
  elif units == "MeV":
    m_mpart = e_per_mparte*m_e
  else:
    raise Exception("The desired units for the momenta are not supported.")

  #Sample the energy of the emitted electrons and get their speed.
  rand_n = np.empty(n)
  for i in range(n):
    rand_n[i] = np.random.uniform()
    while rand_n[i] == 0: #Prevent problems since we need to have some energy beyond E_fermi and E_work.
      rand_n[i] = np.random.uniform()
  E_emitted_electron = np.sqrt(rand_n)*(E_photon-E_work)+E_fermi+E_work
  v = np.sqrt(2*E_emitted_electron/m_e)*c #maginitude of non-relativistic velocity

  #Sample the azimuthal direction of the emitted electrons.
  theta_max = np.arccos(np.sqrt((E_fermi+E_work)/E_emitted_electron))
  theta = np.arccos(1. - (1. - np.sqrt((E_fermi+E_work)/E_emitted_electron))*np.random.uniform(size=n))

  #Sample the polar direction of the emitted electrons
  phi = np.random.uniform(size=n,high=2*np.pi)

  #Calculate the momenta.
  px = m_mpart*v*np.sin(theta)*np.cos(phi)
  py = m_mpart*v*np.sin(theta)*np.sin(phi)
  bz = v/c*np.cos(theta)
  bz2 = np.sqrt(bz**2 - 2.*(E_fermi+E_work)/m_e) #Subtract off the work
  pz = m_mpart*c*bz2

  return  px, py, pz
  
if __name__ == "__main__":
  #Handle arguments
  parser = argparse.ArgumentParser(description='Generates electron macroparticle phase coordinates according to the three setp model.  Output is t,x,y,px,py,pz.')
  parser.add_argument('n', type=float, help='The number of macroparticles desired.')
  parser.add_argument('-n','--number_of_electrons_per_macroparticle', dest="e_per_mpart", type=int, help='The number of electrons simulated by each macroparticle. Default is 100.', default=100)
  parser.add_argument('-t','--std-t', dest="std_t", type=float, help='The standard deviation of the gaussian of electron emission in time. Default is about 21 fs.', default=50./( 2. * np.sqrt( 2. * np.log(2.) ) )*10**(-15))
  parser.add_argument('-R', '-r', dest="R", type=float, help='The radial length.  For the uniform distribution, this is the radius.  For the Guassian, the standard deviation.  Default is about 96 um.', default=np.sqrt(81*115)*1e-6)
  parser.add_argument('-p','--emission_profile', dest="profile", type=str, help='Describes the profile of the electrons in the radial direction.  Options are Gaussian or uniform.  Default is Gaussian.', default="Gaussian")
  parser.add_argument('-o','--output_file', dest="output_file", type=str, help='The path to where the output will be written.  Defualt is IniConditions.txt in the current working directory.', default="IniConditions.txt")
  parser.add_argument('--fermi_energy', dest="E_fermi", type=float, help='The fermi energy of the electron source. Default is 5.000E-6.', default=5.000E-6)
  parser.add_argument('--work_potential', dest="E_work", type=float, help='The work potential of the electron source. Default is 4.450E-6.', default=4.450E-6)
  parser.add_argument('--photon_energy', dest="E_photon", type=float, help='The energy of the lazer\'s photon. Default is 4.650E-6.', default=4.650E-6)
  parser.add_argument('-u','--units', dest="units", type=str, help='Specifies the units of the momenta.  Default is m/s (mks).  Supports mks and MeV/c (specify "MeV").', default="mks")
  parser.add_argument('-d','--delimiter', dest="delimiter", type=str, help='Specifies the delimiter used in the output file.  Default is " ".', default=" ")
  
  args = parser.parse_args()
  args.n = int(args.n)

  #Produce phase space 
  
  t = np.random.normal(scale=args.std_t,size=args.n)
  t = np.sort(t)

  x, y = spatial_transverse_distribution(args.n,args.R,d=args.profile)

  px, py, pz = photoemission_p(**args.__dict__)

  output = np.array([t,x,y,px,py,pz]).transpose()

  with open(args.output_file,'w') as f:
    for i in range(args.n):
      f.write(args.delimiter.join(map(str, output[i])))
      f.write("\n")
