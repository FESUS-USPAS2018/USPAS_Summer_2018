def dump_phase_volume(step,obj,mass):
  """
  Prints out the phase volume of a species to a file named with step.
  Args:
    obj: A container holding the species for which the phase volume will be
      obtained. 
    mass: The mass of the macroparticle used to do momentum conversion.
    step: The iteration after which we are running the dump.  The output file
      will be named step-warp_uem.txt in the running directory.
  Return value:
    None --- but prints to the output file 
  """
  n = obj.getn()
  x = obj.getx()
  y = obj.gety()
  z = obj.getz()
  px = mass*obj.getux()
  py = mass*obj.getuy()
  pz = mass*obj.getuz()
  vx = obj.getvx()
  vy = obj.getvy()
  vz = obj.getvz()

  with open(str(step+1)+"-warp_uem.txt","w") as f:
    for i in range(n):
      output = []
      output.append(x[i])
      output.append(y[i])
      output.append(z[i])
      output.append(px[i])
      output.append(py[i])
      output.append(pz[i])
      output.append(vx[i])
      output.append(vy[i])
      output.append(vz[i])
      output = [str(o) for o in output]
      f.write(" ".join(output) + "\n")
  return 

def dump_phase_volume_with_zlimits(step,obj,mass,zmin=0,zmax=1,**kwargs):
  """
  Prints out the phase volume of a species to a file named with step.
  Args:
    obj: A container holding the species for which the phase volume will be
      obtained. 
    mass: The mass of the macroparticle used to do momentum conversion.
    step: The iteration after which we are running the dump.  The output file
      will be named step-warp_uem.txt in the running directory.
  Return value:
    None --- but prints to the output file 
  """
  n = obj.getn()
  x = obj.getx()
  y = obj.gety()
  z = obj.getz()
  px = mass*obj.getux()
  py = mass*obj.getuy()
  pz = mass*obj.getuz()
  vx = obj.getvx()
  vy = obj.getvy()
  vz = obj.getvz()

  with open(str(step+1)+"-warp_uem.txt","w") as f:
    for i in range(n):
      current_z = z[i]
      if current_z < zmin:
        continue
      if current_z > zmax:
        continue
      output = []
      output.append(x[i])
      output.append(y[i])
      output.append(z[i])
      output.append(px[i])
      output.append(py[i])
      output.append(pz[i])
      output.append(vx[i])
      output.append(vy[i])
      output.append(vz[i])
      output = [str(o) for o in output]
      f.write(" ".join(output) + "\n")
  return 
