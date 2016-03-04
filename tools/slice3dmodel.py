#!/usr/bin/env python3
import os
import sys
import math

listout = []
dict3dcellidto1dcellid = {}
outcellid = 0
with open('model.txt', 'r') as fmodelin:
  npts_model3d = fmodelin.readline()
  t_model = fmodelin.readline() #days
  vmax = fmodelin.readline() #cm/s

  while True:
    blockline1 = fmodelin.readline()
    blockline2 = fmodelin.readline()

    if not blockline1 or not blockline2:
      break

    line1split = blockline1.split()

    if len(line1split) == 5:
      (cellid, posz, posy, posx, rho) = line1split
    else:
      print("WTF? Wrong line size")
      sys.exit()

    line2split = blockline2.split()

    if len(line2split) == 5:
      (ffe, fni, fco, f52fe, f48cr) = map(float,line2split)
    else:
      print("WTF? Wrong line size")
      sys.exit()

    if posz == "0.0000000" and posy == "0.0000000" and float(posx) >= 0.:
      outcellid += 1
      dict3dcellidto1dcellid[int(cellid)] = outcellid
      vel = math.sqrt(float(posx) ** 2 + float(posy) ** 2 + float(posz) ** 2)
      vcell = float(vel) / float(t_model) / 86400 / 1e5
      listout.append("{:6d}  {:8.0f}  {:8.5f}  {:.5f}  {:.5f}  {:.5f}  {:.5f}  {:.5f}".format(outcellid,vcell,math.log10(float(rho)),ffe,fni,fco,f52fe,f48cr))
      print("Cell {:4d} input1: {:}".format(outcellid,blockline1.rstrip()))
      print("Cell {:4d} input2: {:}".format(outcellid,blockline2.rstrip()))
      print("Cell {:4d} output: {:}".format(outcellid,listout[-1]))

with open('model-1dslice.txt', 'w') as fmodelout:
  fmodelout.write("{:7d}\n".format(outcellid))
  fmodelout.write(t_model)
  for line in listout:
    fmodelout.write(line + "\n")

with open('abundances.txt', 'r') as fabundancesin, open('abundances-1dslice.txt', 'w') as fabundancesout:
  currentblock = []
  keepcurrentblock = False
  for line in fabundancesin:
    linesplit = line.split()

    if len(currentblock) + len(linesplit) >= 30:
      if keepcurrentblock:
        fabundancesout.write("  ".join(currentblock) + "\n")
      currentblock = []
      keepcurrentblock = False

    if len(currentblock) == 0:
      currentblock = linesplit
      if int(linesplit[0]) in dict3dcellidto1dcellid.keys():
        outcellid = dict3dcellidto1dcellid[int(linesplit[0])]
        currentblock[0] = "{:6d}".format(outcellid)
        keepcurrentblock = True
    else:
      currentblock.append(linesplit)

if keepcurrentblock:
  print("WARNING: unfinished block")