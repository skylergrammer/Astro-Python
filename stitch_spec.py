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


def interpSpec(wave, spectrum, grid=0.25):
  '''
  Interpolate the spectra to the same resolution.
  '''
  new_wave = np.arange(np.ceil(min(wave)), np.floor(max(wave)), grid)
  
  f = interp1d(wave, spectrum)

  new_spec = f(new_wave)

  return (new_wave, new_spec)


def getVertOffset(wave1, spectrum1, wave2, spectrum2):

  bounds1 = [min(wave1), max(wave1)]
  bounds2 = [min(wave2), max(wave2)]

  overlap1 = np.where((wave1 >= bounds2[0]) & (wave1 <= bounds2[1]))
  overlap2 = np.where((wave2 >= bounds1[0]) & (wave2 <= bounds1[1]))

  if overlap1 and overlap2:
  
    med_val1 = np.median(spectrum1[overlap1])
    med_val2 = np.median(spectrum2[overlap2])

    vert_offset = med_val1 - med_val2

  else:
    vert_offset = 0.0

  return vert_offset
  

def stitchSpectra(wave1, spectrum1, wave2, spectrum2, grid=0.25):
  '''
  Stitch together the two spectra.  Overlapping regions are mean combined.
  '''
  lo_wave = min([wave1[0], wave2[0]])
  hi_wave = max([wave1[-1], wave2[-1]])

  all_wave = np.arange(lo_wave, hi_wave, grid)

  waves = np.array(list(wave1) + list(wave2))
  spec = np.array(list(spectrum1) + list(spectrum2))
 
  stitch_spec = []

  for each in all_wave:
    x = np.where(waves == each)
    # If there are matching wavelengths, fill with mean value
    if x[0].any():
      stitch_spec.append(np.mean(spec[x]))
    # If not, just give it a 1.0
    else:
      stitch_spec.append(0.0)

  return (all_wave, np.array(stitch_spec))
  

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
  wave2, raw_spectrum2, header2 = readSpectrum(filename2)
  smooth_spectrum2 = smoothSpec(raw_spectrum2)

  # Determine vertical offset
  vert_offset = getVertOffset(wave1, smooth_spectrum1, wave2, smooth_spectrum2)

  # Apply vertical offset
  spectrum2 = raw_spectrum2 + vert_offset

  # Verify that the two spectra are the same object
  if header1['object'] != header2['object']:
    sys.exit('ERROR: Input spectra are not the same object!')

  # Determine the interpolation resolution: uses the highest resolution between the spectra
  interp_bin = np.round(min([header1['cdelt1'], header1['cdelt1']]),1)
 
  # Interpolate both spectra to the same resolution
  interp_wave1, interp_spectrum1 = interpSpec(wave1, spectrum1, grid=interp_bin)
  interp_wave2, interp_spectrum2 = interpSpec(wave2, spectrum2, grid=interp_bin)

  # Stitch together the spectra
  x,y = stitchSpectra(interp_wave1, interp_spectrum1, 
                      interp_wave2, interp_spectrum2, grid=interp_bin) 

  # Write stitched spectrum to file
  writeSpectrum(x,y, header1, grid=interp_bin)



def main():

  parser = argparse.ArgumentParser(description='Stitches together two normalized spectra.')
  parser.add_argument('--blue', help='Filename for blue spectrum.')
  parser.add_argument('--red', help='Filename for red spectrum.')
  parser.add_argument('--list', help='Two-column list of spectra to be combined.  Each row must correspond to a single object')
  args = parser.parse_args()

  if args.list:
    f = open(args.list, 'r')
  
    for line in f:
      blue, red = line.split()     
      run(blue, red)

  else:
    run(args.blue, args.red)
  

if __name__ == '__main__':
  main()
