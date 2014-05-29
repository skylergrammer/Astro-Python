#!/usr/bin/env python

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import argparse


def get_ref_data(filename):
  """Reads the file that contains the reference information for each object -- e.g. id, light curve filenames, reference mags, etc.
  """
  header = ['x_shifted','y_shifted','id',   #the header column names
            'lcR','xR','yR','magR','emagR',
            'lcV','xV','yV','magV','emagV',
            'lcB','xB','yB','magB','emagB',
            'lcU','xU','yU','magU','emagU']

  fmt = ['f8','f8','S30', #format keys for each column
         'S40','f8','f8','f8','f8',
         'S40','f8','f8','f8','f8',
         'S40','f8','f8','f8','f8',
         'S40','f8','f8','f8','f8',]

  key = np.genfromtxt(filename, names=header, autostrip=True, dtype=fmt)

  return key


def gen_lc(ref_data, sourceid, zps, dm, show="screen"):
  """Given the reference data, will find the sourceid and locate all the files containing the UBVR light curves.  Feeds the light curve to plot_lc to do the actual plotting.
  """
  #If distance modulus not set, set it to that of M101
  if not dm: dm = 29.05
  
  #Check to see if provided ID is actually in the reference file
  all_ids = [x['id'] for x in ref_data]
  if not sourceid in all_ids:
    sys.exit('%s not in reference file' %sourceid)

  #Get the reference data corresponding to provided ID
  ref_line = [x for x in ref_data if sourceid in x]
  ref_line = ref_line[0]

  if os.path.isfile(ref_line['lcB']) and os.path.isfile(ref_line['lcV']): 
    header = ['mjd', 'dCounts', 'edCounts']
    lcB = np.genfromtxt(ref_line['lcB'], usecols=[0,1,2], names=header)
    lcB_m = convertMags((ref_line['magB'], ref_line['emagB']), sourceid, lcB, dm, 
                       (zps['B'],zps['eB']))
    
    lcV = np.genfromtxt(ref_line['lcV'], usecols=[0,1,2], names=header)
    lcV_m = convertMags((ref_line['magV'], ref_line['emagV']), sourceid, lcV, dm, 
                       (zps['V'],zps['eV']))

    plotTCMD(lcB_m, lcV_m)
  

def convertMags(ref_info, sourceid, lc, dm, zp_info):
  
  ref_mag, ref_mag_err = ref_info
  zp, ezp = zp_info
  offset = zp + 25.0
  
  #convert ref mag to ref counts
  ref_counts = 10**((ref_mag-offset + dm)/-2.5)
  
  #convert dcounts to a magnitude
  counts_i = [ref_counts - x if x < ref_counts else np.nan 
              for x in lc['dCounts']]
  mag_i = [-2.5*np.log10(x)+offset for x in counts_i]
  Beta = (2.5*ref_mag_err)**2
  emag_i = [np.sqrt(Beta + (1.1*y/x)**2 + ezp**2) 
            for x,y in zip(counts_i,lc['edCounts'])]

  return (np.array(lc['mjd']), np.array(mag_i), np.array(emag_i))


def plotTCMD(lc1, lc2):

  fig, ax = plt.subplots(figsize=(14,8), dpi=72)
  fig.subplots_adjust(wspace=0, left=0.125, right=0.98,
                        bottom=0.125, top=0.925)
  t1, m1, em1 = lc1
  t2, m2, em2 = lc2

  exit()
  for each in t1: print each, np.round(each, 1)
  print '           '  
  for each in t2: print each, np.round(each, 1)
  exit()
  plt.plot(m1-m2, m2, 'k')
  plt.ylim(max(m2), min(m2))
  plt.show()

  
 
  
def main():

  parser = argparse.ArgumentParser()
  parser.add_argument('ref', metavar='reference_file')
  parser.add_argument('--chip', choices=['1','2','3','4'], required=True, help='The chip number.  This is important for determining the photometric zeropoints.')
  parser.add_argument('--id', type=str,  help='Options: single source id or "list". "list" will list all the sources, then quit".')
  parser.add_argument('--dm', type=float, help='Distance modulus.  Default is set to 29.05 (M101).')  
  parser.add_argument('--show', choices=["screen", "pdf", "png", "jpg"], default="screen" , help='How the light curve will be viewed.')

  args = parser.parse_args()

  #The zeropoints for each filter
  if args.chip == '1':  
    zps = {'U':7.02, 'eU':0.01, 
           'B':8.65, 'eB':0.01,
           'V':9.22, 'eV':0.01,
           'R':8.74, 'eR':0.01}
  if args.chip == '2':  
    zps = {'U':6.95, 'eU':0.01, 
           'B':9.43, 'eB':0.07,
           'V':9.42, 'eV':0.01,
           'R':8.94, 'eR':0.03}
  if args.chip == '3':  
    zps = {'U':7.00, 'eU':0.01, 
           'B':8.62, 'eB':0.01,
           'V':9.11, 'eV':0.01,
           'R':8.74, 'eR':0.01}
  if args.chip == '4':  
    zps = {'U':7.02, 'eU':0.01, 
           'B':8.67, 'eB':0.01,
           'V':9.17, 'eV':0.01,
           'R':8.78, 'eR':0.01}

  #Retrieves the reference data
  ref_data =  get_ref_data(args.ref)  #read the ref file

  #Lists all the source IDs in the reference data
  if args.id.lower() == "list": 
    for line in ref_data['id']: print line
    exit()
  else:
    gen_lc(ref_data, args.id, zps, args.dm, show=args.show)

if __name__ == '__main__':
  main()
