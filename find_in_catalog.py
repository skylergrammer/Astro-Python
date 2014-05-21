#! /usr/bin/env python
#mod-date:2014-04-18
#upload-date:

import numpy as np
from numpy.lib import recfunctions as rfn
from scipy.constants import *
import matplotlib.pyplot as plt
import pyfits
import argparse
import sys
from sex2deg import sex2deg

def get_targs(targ_filename):
  '''  
  Opens the target file using numpy genfromtxt.  Assumes the format is given by thehead.  If coordinates are sexagesimal then converts them to degrees and appends a two new columns ('degra' and 'degdec') to the record array.  Even if coordinates are already in degrees, the 'degra' and 'degdec' columns are appended but will be identical to the original 'ra' and 'dec'.
  '''  


  f = open(targ_filename,'r')
  header = f.readline().lower().split()
  f.close()

  thehead = ['ra', 'dec', 'lbtid']
  dtype = ['S30', 'S30', 'S30']
 
  if header == thehead:
    data = np.genfromtxt(targ_filename, names=header, dtype=dtype, skip_header=True)
    
    # If the coordinates are in sexagesimal, converts them to degrees
    if ':' in data['ra'][0]:
      ra_degrees = np.array([sex2deg(x) for x in data['ra']])
      dec_degrees = np.array([sex2deg(x, RA=False) for x in data['dec']])
     
      targs = rfn.append_fields(data, names=['degra','degdec'], data=[ra_degrees, dec_degrees], 
                                dtypes=['f8','f8'], usemask=False)

    else:
      targs = rfn.append_fields(data, names=['degra','degdec'], data=[data['ra'], data['dec']], 
                                dtypes=['f8','f8'], usemask=False)
          
    return targs
  
  else:
    sys.exit('Columns need to be title: %s' % thehead)


def get_cat(cat_filename, cut=0):
  '''  
  Opens the target file using numpy genfromtxt.  Assumes the format is given by thehead.  If coordinates are sexagesimal then converts them to degrees and appends a two new columns ('degra' and 'degdec') to the record array.  Even if coordinates are already in degrees, the 'degra' and 'degdec' columns are appended but will be identical to the original 'ra' and 'dec'.
  '''

  f = open(cat_filename,'r')
  header = f.readline().lower().split()
  f.close()

  thehead = ['hstid', 'field', 'ra', 'dec', 'v', 'verr', 'bvcol', 'bvcolerr', 'vicol', 'vicolerr']
  dtype = ['S30', 'S30', 'S30', 'S30', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8']
  
  if header == thehead:
    data = np.genfromtxt(cat_filename, names=header, dtype=dtype, skip_header=True, skip_footer=cut)
  
    if ':' in data['ra'][0]:
      ra_degrees= np.array([sex2deg(x) for x in data['ra']])
      dec_degrees = np.array([sex2deg(x, RA=False) for x in data['dec']])      
      
      cat = rfn.append_fields(data, names=['degra','degdec'], data=[ra_degrees, dec_degrees], 
                              dtypes=['f8','f8'], usemask=False)
    
    else:
      cat = rfn.append_fields(data, names=['degra','degdec'], data=[data['ra'], data['dec']], 
                              dtypes=['f8','f8'], usemask=False)
    
    return cat
  
  else:
    sys.exit('Columns need to be title: %s' % thehead)


def find_nearby(targets, catalog, max_dist=0.5):
  '''
  Finds all the sources within a max distance of each target.  
  '''

  # Determines the squared separation between targets
  sep = lambda x0,y0,x,y: np.power((x-x0)*np.cos(y0*pi/180.),2) + np.power((y-y0),2)
  
  # Convert max distance from arcseconds to squarded degrees
  max_dist = max_dist / 3600.**2
  

  f = open('FIN_output.cat', 'w+')
  ff = open('FIN_faintsources.cat', 'w+')
  fff = open('FIN_notmatched.cat', 'w+')

  for each in targets:

    rvec = sep(each['degra'], each['degdec'], catalog['degra'], catalog['degdec'])

    nearby = np.argsort(rvec)

    if any(nearby) and rvec[nearby[0]] <= max_dist: 
      jj = catalog[nearby[0]]
      fill = [jj['hstid'], each['lbtid'], jj['ra'], jj['dec'], jj['v'], jj['bvcol'], jj['vicol'], '\n']
      f.write('   '.join([str(x) for x in fill]))
      if jj['v'] > 23.5:
        ff.write('   '.join([str(x) for x in fill]))
    else:
      fff.write('   '.join([str(x) for x in each])+'\n') 
      fill0 = ['None', each['lbtid'], each['ra'], each['dec'], 'None', 'None', 'None', '\n']
      f.write('   '.join([str(x) for x in fill0]))

  f.close()
  ff.close()
  fff.close()

def main():
  
  parser = argparse.ArgumentParser()
  parser.add_argument('--t', metavar='targlist', required=True)
  parser.add_argument('--c', metavar='catalog', required=True)
  parser.add_argument('--tol', metavar='tolerance', type=float, default=0.5, help='Max separation in arcseconds. Default is 0.5".')

  args = parser.parse_args()

  targs = get_targs(args.t)
  cat = get_cat(args.c, cut=512099)
  
  find_nearby(targs, cat, max_dist=args.tol)


if __name__ == '__main__':
  main()
