import numpy as np
import cPickle as pickle
from warp import *

def phase_volume_pickle_loader(pickle_dict_file,time_conversion=1.,
          position_conversion=1.,momentum_conversion=1.,**kwargs):
  """
  Read in the initial conditions data that is stored in a pickled
  dict.
  Args:
    pickle_dict_file: The pickle file from which we will read.
    *_conversion: Will multiply the corresponding coordinates
      to convert them to appropriate coordinates
  Return value:
    (t,x,y,z,px,py,pz) : The phase coordinates for the N particles
      stored as single dimensioned np arrays. 
  """
  ff = open(pickle_dict_file, 'r')
  dd = pickle.load(ff)
  ff.close() 

  # --- create ... arrays to store data  
  n  = len(dd)
  t  = np.zeros(n)  
  x  = np.zeros(n) 
  y  = np.zeros(n) 
  z  = np.zeros(n) 
  px = np.zeros(n) 
  py = np.zeros(n) 
  pz = np.zeros(n) 

  # --- place data in ... arrays.   
  for ii in range(n):
    row = dd[ii]
    t[ii] = row['t']*time_conversion
    #
    x[ii] = row['x']*position_conversion
    y[ii] = row['y']*position_conversion
    z[ii] = row['z']*position_conversion
    # 
    px[ii] = row['px']*momentum_conversion 
    py[ii] = row['py']*momentum_conversion
    pz[ii] = row['pz']*momentum_conversion

  return (t,x,y,z,px,py,pz)

def photoemission_loader(filepath,time_conversion=1.,
          position_conversion=1.,momentum_conversion=1.,
          delimiter=" ",**kwargs):
  """
  Read in the initial conditions data that is stored in a 
  file with the columns t,x,y,px,py,pz
  Args:
    filepath: The file from which we will read.
    *_conversion: Will multiply the corresponding coordinates
      to convert them to appropriate coordinates
  Return value:
    (t,x,y,z,px,py,pz) : The phase coordinates for the N particles
      stored as single dimensioned np arrays, with z set to zeros. 
  """
  t,x,y,px,py,pz = np.loadtxt(filepath,delimiter=delimiter).transpose()

  z  = np.zeros(t.size) 

  #Do conversions
  t = time_conversion*t
  x = position_conversion*x
  y = position_conversion*y
  px = momentum_conversion*px
  py = momentum_conversion*py
  pz = momentum_conversion*pz

  return (t,x,y,z,px,py,pz)

def get_pulse_velocity_from_momentum(coordinate_array_dict,mass,direction="z"):
  """
  Calculates the boost to the mean "particle" assuming the mean
  momentums in the other two directions are 0. 
  Args:
    coordinate_array_dict: A dictionary with at least the keys
      z, px, py, pz with values for the corresponding np arrays. 
    mass: Mass of the particle to rescale the momentum.
    direction: The desired direction for which the velocity will be calculated.
  Return value:
    The mean velocity in the provided direction.
  """
  clight = 299792458
  mean_p = np.mean(coordinate_array_dict["p"+direction])
  gamma  = np.sqrt(1. + (mean_p/(mass*clight))**2 )
  return mean_p/(mass*gamma)

def get_distribution_average(filepath,component="z"):
  """
  Reads in the distbution coordinates and calculates the returns the mean
  z value.
  Args:
    filepath:  The path to the particles coordinates.
    component: The desired over which the average will be taken.
  Return value:
    The mean position in the provided direction.
  """
  try:
    keys = ["x", "y", "z", "3", "4", "5", "6", "7", "8"]
    data = getdatafromtextfile(filepath,nskip=0,dims=[9,None]) 
  except:
    keys = ["x", "y", "z", "3", "4", "5"]
    data = getdatafromtextfile(filepath,nskip=0,dims=[6,None]) 
  return data[keys.index(component)].mean()
