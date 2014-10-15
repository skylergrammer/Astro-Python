#!/usr/bin/env python

import argparse
from scipy.ndimage.filters import *
import numpy as np
import pyfits
import matplotlib.pyplot as plt

def read_spectrum(filename, units='pixels'):
  '''
  Read in the spectrum and output the x values, as pixels or angstroms, and the y
  values as whatever are read in.
  '''
  data, header = pyfits.getdata(filename, header=True)
  
  if units == 'pixels':
    x = np.arange(0, header['naxis1'])

  else:
    lambda1 = header['crval1']
    lambda2 = float(header['naxis1'])*header['cdelt1']+float(header['crval1'])
  
    # Convert pixels to angstroms
    x = np.arange(lambda1, lambda2, header['cdelt1'])
    
  return (x, data[0], header)


def smooth(data, header, size, units='pixels'):

  if units == 'physical':
    size = size * header['cdelt1']
    
  smoothed_data = median_filter(data, size=size)

  return smoothed_data


def write_spectrum(data, header, orig_fn):

  new_fn = orig_fn.replace('.fits', '_smooth.fits')

  pyfits.writeto(new_fn, data, header=header, clobber=True)
 
  return new_fn 


def main():

  parser = argparse.ArgumentParser()
  parser.add_argument('spectra', nargs='+')
  parser.add_argument('--size', type=float, default=3)
  parser.add_argument('--units', choices=['physical','pixels'], default='pixels')

  args = parser.parse_args()

  for each in args.spectra:

    x, y, header = read_spectrum(each, units=args.units)

    y_smooth = smooth(y, header, args.size, units=args.units)

    smoothed_spec_fn = write_spectrum(y_smooth, header, each)

if __name__ == "__main__":
  main()