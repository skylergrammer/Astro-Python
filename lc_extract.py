#! /usr/bin/env python
import sys
import os
import numpy as np
import argparse
import matplotlib.pyplot as plt

def get_info_table(filename):
  '''
  Reads the file that contains the reference information for each object 
  -- e.g. id, light curve filenames, reference mags, etc.
  '''

  # Header column names
  header = ['x_shifted','y_shifted','id',   
            'lcR','xR','yR','magR','emagR',
            'lcV','xV','yV','magV','emagV',
            'lcB','xB','yB','magB','emagB',
            'lcU','xU','yU','magU','emagU']

  # Format for each column
  fmt = ['f8','f8','S30', 
         'S40','f8','f8','f8','f8',
         'S40','f8','f8','f8','f8',
         'S40','f8','f8','f8','f8',
         'S40','f8','f8','f8','f8',]

  info_table = np.genfromtxt(filename, names=header, autostrip=True, dtype=fmt,
                             missing_values=("0.0","000"), filling_values=np.nan)

  return info_table


class ZeroPoints: 
  def __init__(self, data, dm):
    self.U = data["U"]
    self.eU = data["eU"]
    self.B = data["B"]
    self.eB = data["eB"]
    self.V = data["V"]
    self.eV = data["eV"]
    self.R = data["R"]
    self.eR = data["eR"]
    self.dm = dm[0]
    self.edm = dm[1]
    
    
def get_zeropoints(filename):
  '''
  Reads the zero point info from the zero point file.
  '''

  header = ["chip","U","eU","B","eB","V","eV","R","eR"]
  fmt = ['L', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8']
 
  zeropoints = np.genfromtxt(filename, names=header)

  return zeropoints


def get_lc_data(filename):
  '''
  Reads in the actual light curve info from the light curve filename provided
  in the info table.
  '''

  header = ['mjd', 'dCounts', 'edCounts']
  lc_data = np.genfromtxt(filename, usecols=[0,1,2], names=header)

  return lc_data


class RefInfo:
  '''
  Class to hold all the reference info.  Includes reference mags for each 
  bandpass, errors, and the filename for the light curve.
  '''

  def __init__(self, data, id):
    self.id = id
      
    ref_info = [x for x in data if id in x]

    if not ref_info:
      sys.exit("\n***ERROR: %s does not exist within the info table.\n" % id)
    else:
      ref_info = ref_info[0]
    
    if os.path.isfile(ref_info["lcU"]):
      self.lcU = ref_info["lcU"]
      self.refU = ref_info["magU"]
      self.erefU = ref_info["emagU"]
    else:
      self.lcU = False
      self.refU = False
      self.erefU = False

    if os.path.isfile(ref_info["lcB"]):
      self.lcB = ref_info["lcB"]
      self.refB = ref_info["magB"]
      self.erefB = ref_info["emagB"]
    else:
      self.lcB = False
      self.refB = False
      self.erefB = False

    if os.path.isfile(ref_info["lcV"]):
      self.lcV = ref_info["lcV"]
      self.refV = ref_info["magV"]
      self.erefV = ref_info["emagV"]
    else:
      self.lcV = False
      self.refV = False
      self.erefV = False

    if os.path.isfile(ref_info["lcR"]):
      self.lcR = ref_info["lcR"]
      self.refR = ref_info["magR"]
      self.erefR = ref_info["emagR"]
    else:
      self.lcR = False
      self.refR = False
      self.erefR = False

  def lc_exist(self):
    '''
    Call to print the file name for each light curve.  If it exists, prints a 
    file name.  If it does not exist, then indicates no light curve exists.
    '''

    if self.lcU: 
      print(self.lcU)
    else: 
      print("No U light curve")
    if self.lcB:
      print(self.lcB)
    else: 
      print("No B light curve")
    if self.lcV:
      print(self.lcV)
    else: 
      print("No V light curve")
    if self.lcR:
      print(self.lcR)
    else: 
      print("No R light curve")


def convert_to_lc(lc, src_ref_info, zps, band=None, chip=1):
  '''
  This function will convert counts, with respect to reference magnitude, to 
  magnitudes for each epoch and determine the error.
  '''

  # chips values range 1:4 but chip indices 0:3
  chip = chip - 1

  # Set the reference mags based on the given bandpass  
  if band == "U":
    ref_mag = src_ref_info.refU
    e_ref_mag = src_ref_info.erefU
    zp = zps.U[chip]
    ezp = zps.eU[chip]

  elif band == "B":
    ref_mag = src_ref_info.refB
    e_ref_mag = src_ref_info.erefB
    zp = zps.B[chip]
    ezp = zps.eB[chip] 

  elif band == "V":
    ref_mag = src_ref_info.refV
    e_ref_mag = src_ref_info.erefV
    zp = zps.V[chip]
    ezp = zps.eV[chip] 

  elif band == "R":
    ref_mag = src_ref_info.refR
    e_ref_mag = src_ref_info.erefR
    zp = zps.R[chip]
    ezp = zps.eR[chip] 

  # Default offset is the zeropoint plus the IRAF 25.0
  offset = zp + 25.0

  
  # Convert ref_mag to counts
  ref_counts = 10**((ref_mag-offset + zps.dm)/-2.5)
  
  # Convert delta counts to a magnitude
  counts_i = [ref_counts - x if x < ref_counts else np.nan 
              for x in lc['dCounts']]
  mag_i = np.array([-2.5*np.log10(x)+offset for x in counts_i])
  
  # Factor to be used in error propagation
  Beta = (2.5*e_ref_mag)**2
  
  # Errors in magnitudes
  emag_i = [np.sqrt(Beta + (1.1*y/x)**2 + ezp**2) 
            for x,y in zip(counts_i,lc['edCounts'])]  
  
  # Output dictionary to hold light curve
  output = {"mjd":lc["mjd"], "mag":mag_i, "emag":emag_i}
  
  return output


def extract_light_curves(src_ref_info, zps, chip=None):
  '''
  Extracts the light curves for each band.  
  '''

  if src_ref_info.lcU:
    lcU_data = get_lc_data(src_ref_info.lcU)
    light_curve_U = convert_to_lc(lcU_data, src_ref_info, zps, band="U", chip=chip)
  

  if src_ref_info.lcB:
    lcB_data = get_lc_data(src_ref_info.lcB)
    light_curve_B = convert_to_lc(lcB_data, src_ref_info, zps, band="B", chip=chip)

  if src_ref_info.lcV:
    lcV_data = get_lc_data(src_ref_info.lcV)
    light_curve_V = convert_to_lc(lcV_data, src_ref_info, zps, band="V", chip=chip)

  if src_ref_info.lcR:
    lcR_data = get_lc_data(src_ref_info.lcR)
    light_curve_R = convert_to_lc(lcR_data, src_ref_info, zps, band="R", chip=chip)


def main():
  
  parser = argparse.ArgumentParser()
  parser.add_argument("--info", required=True)
  parser.add_argument("--source", required=True)
  parser.add_argument("--zp", required=True)
  parser.add_argument("--dm", nargs=2, type=float, required=True)
  parser.add_argument("--chip", type=int, required=True)
  args = parser.parse_args()
  
  info_table = get_info_table(args.info)
  
  raw_zp = get_zeropoints(args.zp)
  zps = ZeroPoints(raw_zp, args.dm)

  src_ref_info = RefInfo(info_table, args.source)    
  src_ref_info.lc_exist()

  extract_light_curves(src_ref_info, zps, chip=args.chip)


if __name__ == "__main__":
  main()