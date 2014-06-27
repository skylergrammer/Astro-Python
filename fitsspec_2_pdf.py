#! /usr/bin/env python

import pyfits
import numpy as np
import matplotlib.pyplot as plt
import argparse
from scipy.ndimage.filters import *
from matplotlib.ticker import MultipleLocator, FormatStrFormatter


def getSpectrum(filename):

  # Retrieve data
  data, header = pyfits.getdata(filename, header=True)

  if len(data) != int(header['naxis1']): data = data[0]

  # Convert pixels to wavelengths
  lambda1 = header['crval1']
  lambda2 = float(header['naxis1'])*header['cdelt1']+float(header['crval1'])
  
  wave = np.arange(lambda1, lambda2, header['cdelt1'])
  wave = wave[:len(data)]
  
  return (wave, data, header)

def smoov(data, size):

  smooved_data = gaussian_filter1d(data, size)

  return smooved_data


def plotSpec(wave, data, header, bounds):

  median_y = np.median(data[np.where((wave > bounds[0]) & (wave <= bounds[1]))])
  max_y = max(data[np.where((wave > bounds[0]) & (wave <= bounds[1]))])
  min_y = min(data[np.where((wave > bounds[0]) & (wave <= bounds[1]))])
  std_y = np.std(data[np.where((wave > bounds[0]) & (wave <= bounds[1]))])


  fig, ax = plt.subplots(figsize=(15,5), dpi=72)
  fig.subplots_adjust(wspace=0.2, left=0.025, right=0.975,
                        bottom=0.2, top=0.92)
  

  plt.plot(wave, data, 'k', lw=1.5, label=header['object'])
  
  # x and y limits
  plt.ylim(median_y - std_y, median_y + std_y)
  plt.xlim(bounds[0], bounds[1])
  
  # x-tick major and minor locations at 500 and 100 Angstroms
  #xmajorLocator = MultipleLocator(500)
  #xminorLocator = MultipleLocator(100)
  #ax.xaxis.set_major_locator(xmajorLocator)
  #ax.xaxis.set_minor_locator(xminorLocator)
  
  # Suppress ticks on y axiss
  plt.tick_params(axis='x', which='both', labelsize=22)
  ax.tick_params(axis='y', which='both', left='off',right='off', labelleft='off')

  # Set x label
  plt.xlabel('$\lambda$ [$\AA$]', fontsize=30)
  ax.xaxis.labelpad = 10
  
  # Set title
  plt.title(header['object'], fontsize=30)
  plt.show()

def main():

  parser = argparse.ArgumentParser()
  parser.add_argument('input', nargs='+')
  parser.add_argument('--bounds', nargs=2, type=float, required=True)
  parser.add_argument('--smooth', type=int, default=2)

  args = parser.parse_args()


  for each in args.input:

    wave, data, header = getSpectrum(each)
    smooved_data = smoov(data, args.smooth)

    plotSpec(wave, smooved_data, header, args.bounds)

if __name__ == '__main__':
  main()
  