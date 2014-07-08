#! /usr/bin/env python

import pyfits
import numpy as np
import matplotlib.pyplot as plt
import argparse
from scipy.ndimage.filters import *
from matplotlib.ticker import MultipleLocator
import sys

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


def plotSpec(wave, data, header, bounds, fig, ax, number):

  median_y = np.median(data[np.where((wave > bounds[0]) & (wave <= bounds[1]))])
  max_y = max(data[np.where((wave > bounds[0]) & (wave <= bounds[1]))])
  min_y = min(data[np.where((wave > bounds[0]) & (wave <= bounds[1]))])
  std_y = np.std(data[np.where((wave > bounds[0]) & (wave <= bounds[1]))])

  # Plot the spectrum
  ax[number].plot(wave, data, 'k', lw=1, label=header['object'])
  
  # x and y limits
  med_blue = np.median(data[0:100])
  med_red = np.median(data[-100:-1])
  
  if med_blue >= med_red:
    ax[number].set_ylim(min_y, 2*med_blue)
  else:
    ax[number].set_ylim(min_y, 2*med_red)

  ax[number].set_xlim(bounds[0], bounds[1])
  

def main():

  parser = argparse.ArgumentParser()
  parser.add_argument('input', nargs='+')
  parser.add_argument('--l', action='store_true')
  parser.add_argument('--bounds', metavar=('Blue','Red'), nargs=2, type=float, required=True)
  parser.add_argument('--smooth', type=int, default=2)
  parser.add_argument('--ydim', metavar='Rows', type=int, default=1)

  args = parser.parse_args()

  xmajorLocator = MultipleLocator(500)
  xminorLocator = MultipleLocator(100)

  if args.l:
    file_list = open(args.input[0], 'r')
  else:
    file_list = args.input
    
  if len(file_list) < 2:
    sys.exit('Must have at least 2 spectra to plot.')      
  
  ax_count = 0
  for each in file_list:

    # If restarting: initialize plot    
    if ax_count == 0:
      fig, ax = plt.subplots(nrows=args.ydim, ncols=1, sharex=True, sharey=False,
                             figsize=(15,2*args.ydim), dpi=72)
      fig.subplots_adjust(hspace=0, left=0.025, right=0.975,
                          bottom=0.15, top=0.99)
      
    # Read in spectrum
    wave, data, header = getSpectrum(each)
    # Smooth spectrum
    smooved_data = smoov(data, args.smooth)
    # Plot spectrum
    plotSpec(wave, smooved_data, header, args.bounds, fig, ax, ax_count)
    ax[ax_count].legend(frameon=False, loc='upper left', markerscale=None) 
    ax[ax_count].xaxis.set_major_locator(xmajorLocator)
    ax[ax_count].xaxis.set_minor_locator(xminorLocator)

    # If bottom plot, make the xlabel
    if ax_count == args.ydim-1:
      ax[ax_count].tick_params(axis='y', which='both', left='off',right='off', labelleft='off')
      ax[ax_count].tick_params(axis='x', which='both', labelsize=22)
      ax[ax_count].set_xlabel('$\lambda$ [$\AA$]', fontsize=30)
      ax[ax_count].xaxis.labelpad = 10
    
    # If not the bottom plot, suppress all labels
    else: 
      ax[ax_count].tick_params(axis='both', which='both', 
                               bottom='on', top='off', left='off',right='off', 
                               labelleft='off', labelbottom='off')

    # If ax_count is less than the number of rows, increase value by 1
    if ax_count < args.ydim-1:
      ax_count += 1
    # If ax_count equal to number of rows, reset it and show the final plot
    else: 
      ax_count = 0
      plt.show()

  '''  
  plt.savefig('spec_'+header['object']+'.pdf', 
                  dpi=None, facecolor='w', edgecolor='w',
                  transparent=True, bbox_inches=None, pad_inches=0.1,
                  frameon=None)
  '''
if __name__ == '__main__':
  main()
  
