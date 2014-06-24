#! /usr/bin/env python


import numpy as np
import matplotlib.pyplot as plt
import pyfits
import argparse
from scipy.ndimage.filters import gaussian_filter1d


def readSpectrum(filename):
  
  data, header = pyfits.getdata(filename, header=True)

  lambda1 = header['crval1']
  lambda2 = float(header['naxis1'])*header['cdelt1']+float(header['crval1'])
  
  # Convert pixels to angstroms
  wave = np.arange(lambda1, lambda2, header['cdelt1'])
  wave = wave[:len(data[0])]

  return (wave, data[0], header)


def smoothSpec(spectrum):
  
  smooth_spec = gaussian_filter1d(spectrum, 5)

  return smooth_spec


def findMultiplier(wave1, spectrum1, wave2, spectrum2):

  bounds1 = (min(wave1), max(wave1))
  bounds2 = (min(wave2), max(wave2))

  overlap1 = [y for x,y in zip(wave1, spectrum1) if x >= bounds2[0] and x <= bounds2[1]]
  overlap2 = [y for x,y in zip(wave2, spectrum2) if x >= bounds1[0] and x <= bounds1[1]]

  ratio = np.median(overlap1)/np.median(overlap2)
  ratio = np.median(overlap1)-np.median(overlap2)
  return ratio


def stitch(wave1, spectrum1, wave2, spectrum2):

  return
def main():

  parser = argparse.ArgumentParser()
  parser.add_argument('spec1')
  parser.add_argument('spec2')

  args = parser.parse_args()

  wave1, spectrum1, header1 = readSpectrum(args.spec1)
  smooth_spectrum1 = smoothSpec(spectrum1)

  wave2, spectrum2, header2 = readSpectrum(args.spec2)
  smooth_spectrum2 = smoothSpec(spectrum2)

  ratio = findMultiplier(wave1, smooth_spectrum1, wave2, smooth_spectrum2)

  scaled_spectrum1 = smooth_spectrum1-ratio
  
  plt.plot(wave1, scaled_spectrum1, 'k')
  plt.plot(wave2, smooth_spectrum2, 'r')
  plt.xlim(3500, 8500)
  plt.show()

if __name__ == '__main__':
  main()