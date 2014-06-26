#! /usr/bin/env python


import numpy as np
import matplotlib.pyplot as plt
import pyfits
import argparse
from scipy.ndimage.filters import gaussian_filter1d
from scipy.interpolate import *
import sys

def readSpectrum(filename):
  '''
  Read in a spectrum.  Return wavelengths, intensities, and header.
  '''
  data, header = pyfits.getdata(filename, header=True)
  if len(data) != int(header['naxis1']): data = data[0]

  lambda1 = header['crval1']
  lambda2 = float(header['naxis1'])*header['cdelt1']+float(header['crval1'])
  
  # Convert pixels to angstroms
  wave = np.arange(lambda1, lambda2, header['cdelt1'])
  wave = wave[:len(data)]

  return (wave, data, header)


def smoothSpec(spectrum):
  '''
  Smooth the spectrum so that small oscillations do not affect offset between 
  the overlapping regions.
  '''
  smooth_spec = gaussian_filter1d(spectrum, 5)

  return smooth_spec


def findOffset(wave1, spectrum1, wave2, spectrum2):
  '''
  Find the vertical offset between the overlapping regions in the two spectra.
  '''

  bounds1 = (min(wave1), max(wave1))
  bounds2 = (min(wave2), max(wave2))

  # Find overlapping regions
  overlap1 = [y for x,y in zip(wave1, spectrum1) 
              if x >= bounds2[0] and x <= bounds2[1]]

  overlap2 = [y for x,y in zip(wave2, spectrum2) 
              if x >= bounds1[0] and x <= bounds1[1]]

  # Find the offset between the mean values in the offset regions.
  offset = np.median(overlap1)-np.median(overlap2)
  
  return offset


def interpSpec(wave, spectrum, grid=0.25):
  '''
  Interpolate the spectra to the same resolution.
  '''
  new_wave = np.arange(min(wave), max(wave), grid)

  f = interp1d(wave, spectrum)

  new_spec = f(new_wave)

  return (new_wave, new_spec)


def stitchSpectra(wave1, spectrum1, wave2, spectrum2, grid=0.25):
  '''
  Stitch together the two spectra.  Overlapping regions are mean combined.
  '''
  lo_wave = min([wave1[0], wave2[0]])
  hi_wave = max([wave1[-1], wave2[-1]])

  all_wave = np.arange(lo_wave, hi_wave, grid)

  waves = np.array(list(wave1) + list(wave2))
  spec = np.array(list(spectrum1) + list(spectrum2))

  #stitch_spec = np.array([np.mean(spec[np.where(waves == x)]) for x in all_wave])

  return (waves, spec)
  

def writeSpectrum(wave, spectrum, header, grid=0.25):
  '''
  Write out the stitched spectrum to a fits file.  Keeps the header of one 
  spectrum.  Updates crval1 and cdelt1
  '''
  header.update(CRVAL1=min(wave))
  header.update(CDELT1=grid)
  
  pyfits.writeto('spec_'+header['object']+'.fits', spectrum, header=header, clobber=True)


def run(filename1, filename2):

  # Read in spectrum 1 and then smooth it
  wave1, spectrum1, header1 = readSpectrum(filename1)
  smooth_spectrum1 = smoothSpec(spectrum1)
  
  # Read in spectrum 2 and then smooth it
  wave2, spectrum2, header2 = readSpectrum(filename2)
  smooth_spectrum2 = smoothSpec(spectrum2)

  #plt.plot(wave1, smooth_spectrum1, 'k')
  #plt.plot(wave2, smooth_spectrum2, 'r')
  #plt.xlim(3500, 8500)
  #plt.show()

  # Verify that the two spectra are the same object
  if header1['object'] != header2['object']:
    sys.exit('ERROR: Input spectra are not the same object!')

  # Determine the vertical offset between the overlapping regions
  #offset = findOffset(wave1, smooth_spectrum1, wave2, smooth_spectrum2)

  # Subtract the offset between overlapping region
  #mod_spectrum2 = np.array(spectrum2)+offset
  
  # Determine the interpolation resolution: uses the highest resolution between the spectra
  interp_bin = np.round(min([header1['cdelt1'], header1['cdelt1']]),1)
 
  # Interpolate both spectra to the same resolution
  interp_wave1, interp_spectrum1 = interpSpec(wave1, spectrum1, grid=interp_bin)
  interp_wave2, interp_spectrum2 = interpSpec(wave2, spectrum2, grid=interp_bin)

  # Stitch together the spectra
  x,y = stitchSpectra(interp_wave1, interp_spectrum1, 
                      interp_wave2, interp_spectrum2, grid=interp_bin) 

  # Write stitched spectrum to file
  try:
    writeSpectrum(x,y, header1, grid=interp_bin)
  except:
    sys.exit('Cannot write to file.  Verify that output filename does not already exist.')


def main():

  parser = argparse.ArgumentParser(description='Stitches together two overlapping \
                                   spectra.  Use this program if the overlapping \
                                   region is small compared to the full wavelength \
                                   range.')
  parser.add_argument('--blue', help='Filename for blue spectrum.')
  parser.add_argument('--red', help='Filename for red spectrum.')
  parser.add_argument('--list', help='Two-column list of spectra to be combined. \
                                     Each row must correspond to a single object')
  args = parser.parse_args()

  if args.list:
    f = open(args.list, 'r')
    for line in f:
      blue, red = line.split()
      print 'Combining %s with %s' % (blue, red)
      try:
        run(blue, red)
      except:
        continue
  else:
    run(args.blue, args.red)
  

if __name__ == '__main__':
  main()