#! /usr/bin/env python
#mod-date:2014-04-16
#upload-date: 2014-04-16 16:46

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

  thehead = ['hstid', 'lbtid', 'ra', 'dec', 'v', 'bvcol', 'vicol']
  dtype = ['S30', 'S30', 'S30', 'S30', 'f8', 'f8', 'f8']

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


def find_nearby(targets, catalog, max_dist=10):
  '''
  Finds all the sources within a max distance of each target.  
  '''

  # Determines the squared separation between targets
  sep = lambda x0,y0,x,y: np.power((x-x0)*np.cos(y0*pi/180),2) + np.power((y-y0),2)
  
  # Convert max distance from arcseconds to squarded degrees
  max_dist = max_dist / 3600.**2

  for each in targets:

    rvec = sep(each['degra'], each['degdec'], catalog['degra'], catalog['degdec'])


    # Ensure that target is recovered in catalog   
    nearby = np.where(rvec < max_dist)

    if each['hstid'] in catalog['hstid'][nearby]:
      write_cmd_file(catalog[nearby], each)

    else: 
      continue    
      print 'ERROR: Failed to find %s.' % each['hstid']


def write_cmd_file(near_targ, target):
  '''
  Takes the rec array of sources near target and the rec array of the target and produces a fits table.
  '''

  # Columns to be in the fits table: these data are for the nearby sources
  c1 = pyfits.Column(name='HSTID', format='20A', array=near_targ['hstid'])
  c2 = pyfits.Column(name='RA', format='F', array=near_targ['degra'])
  c3 = pyfits.Column(name='DEC', format='F', array=near_targ['degdec'])
  c4 = pyfits.Column(name='V', format='F', array=near_targ['v'])
  c5 = pyfits.Column(name='VERR', format='F', array=near_targ['verr'])
  c6 = pyfits.Column(name='BV', format='F', array=near_targ['bvcol'])
  c7 = pyfits.Column(name='BVERR', format='F', array=near_targ['bvcolerr'])
  c8 = pyfits.Column(name='VI', format='F', array=near_targ['vicol'])
  c9 = pyfits.Column(name='VIERR', format='F', array=near_targ['vicolerr'])

  # Make table
  table_hdu = pyfits.new_table([c1, c2, c3, c4, c5, c6, c7, c8, c9])


  # Updates header with contains the target's info
  table_hdu.header.update(key='HSTID', value=target['hstid'])
  table_hdu.header.update(key='LBTID', value=target['lbtid'])
  table_hdu.header.update(key='RA', value=str(target['ra']))
  table_hdu.header.update(key='DEC', value=str(target['dec']))     

  # Table data cannot be the Primary HDU, so we make an empty Primary HDU
  phdu = pyfits.PrimaryHDU()

  # Zeroth extension is empty, first extension contains the table  
  hdulist = pyfits.HDUList([phdu, table_hdu])
  hdulist.writeto(target['lbtid']+'.fits')

  
def main():
  
  parser = argparse.ArgumentParser()
  parser.add_argument('--t', required=True, metavar='targlist')
  parser.add_argument('--c', required=True, metavar='catalog')
  parser.add_argument('--r', type=float, default=10, metavar='Maximum radial distance in arcseconds.')
  args = parser.parse_args()

  targs = get_targs(args.t)
  cat = get_cat(args.c)
  
  find_nearby(targs, cat, max_dist=args.r)


if __name__ == '__main__':
  main()
