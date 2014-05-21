#! /usr/bin/env python
#mod date: 2014-05-17
#upload date: dickbutt
#version: 1.0

import numpy as np
import matplotlib.pyplot as plt
import argparse
import pyfits
from scipy.interpolate import *
from scipy.ndimage.filters import *
from scipy.stats import *

def getSpectrum(filename):

  data, header = pyfits.getdata(filename, header=True)

  lambda1 = header['crval1']
  lambda2 = float(header['naxis1'])*header['cdelt1']+float(header['crval1'])
  wave = np.arange(lambda1, lambda2, header['cdelt1'])
  
  wave = wave[:len(data[0])]
  
  return (wave, data[0], header)


def normSpectrum(xcoord, ycoord, niter=3, filter_size=50, trim=(100,100)):
  # Trim the input spectrum to remove spurious data
  xcoord = xcoord[trim[0]:-trim[1]]
  ycoord = ycoord[trim[0]:-trim[1]]

  # Smooth the data using a median filter
  ycoord_smooth = median_filter(ycoord, filter_size)
  
  # Apply a Sobel filter to identify lines and then filter it with a 3A Gaussian
  ycoord_sobel = gaussian_filter(sobel(ycoord), 3) 
  
  # The bins used to grid the data; will have bins that are 10% the filter size
  step = 0.1*filter_size
  xi = list(np.arange(min(xcoord), max(xcoord), step))
  
  # Grid the sobel filtered spectrum 
  yi = griddata(xcoord, ycoord_sobel, xi)

  # Now we iterate over the Sobel filtered spectrum to remove lines
  for iteration in range(niter):
    print 'PERFORMING ITERATION %s' %iteration
    
    dyi = yi 
    dyi_smooth = gaussian_filter1d(yi, filter_size)
    f0 = interp1d(xi, dyi_smooth)

    frac_diff = (dyi - f0(xi))/f0(xi)
    mlim = np.std(frac_diff)
    frac_diff_med = tmean(frac_diff, limits=(-3*mlim, 3*mlim))
    frac_diff_std = tstd(frac_diff, limits=(-3*mlim, 3*mlim))
    
    xi_clean = [x for x,y,z, in zip(xi,yi,frac_diff) if np.abs(z) < frac_diff_std]
    yi_clean = griddata(xcoord, ycoord_sobel, xi_clean)

    diagnostics = [(x,z) for x,y,z, in zip(xi,yi,frac_diff) if np.abs(z) > frac_diff_std]

    rej,ect = zip(*diagnostics)
    
    '''
    plt.plot(xi, frac_diff, 'k.')
    plt.plot([min(xi),max(xi)],[frac_diff_med,frac_diff_med],'g')
    plt.plot([min(xi),max(xi)],[-frac_diff_std,-frac_diff_std],'r')
    plt.plot([min(xi),max(xi)],[frac_diff_std,frac_diff_std],'r')
    plt.ylim(max(frac_diff), min(frac_diff))  
    plt.plot(rej, ect, 'r.')    
    plt.show()
    '''

    xi = xi_clean
    yi = yi_clean

  yi = griddata(xcoord, ycoord_smooth, xi_clean)
 
  ind_crop = np.where((xcoord <= max(xi_clean)) & (xcoord >= min(xi_clean)))  
  xcoord_crop = xcoord[ind_crop]
  ycoord_crop = ycoord[ind_crop]
 
  
  f1= interp1d(xi, yi, copy=True, bounds_error=True, fill_value=np.nan)
  yn = f1(xcoord_crop)  
  
  plt.plot(xcoord, ycoord, 'k')
  plt.plot(xcoord, ycoord_smooth, 'g')
  plt.plot(xcoord_crop, yn, 'r')
  plt.show()

  norm_spec = gaussian_filter1d(ycoord_crop/yn, 2)
  plt.plot(xcoord_crop, norm_spec, 'k')
  plt.show()
  return (xcoord_crop, norm_spec)
  
def main():

  parser = argparse.ArgumentParser()
  parser.add_argument('files', nargs='+')
  parser.add_argument('--niter', metavar='niter', type=int, default=3, help='Number of iterations to perform.')
  parser.add_argument('--fs', metavar='Size', type=int, default=100, help='Median filter size to use')
  parser.add_argument('--trim', metavar=('Left','Right'), nargs=2, type=int, default=(100,100), help='Amount to be trimmed off from the edges') 
  args = parser.parse_args()

  for each in args.files:
    wave, data, header = getSpectrum(each)

    wave_norm, data_norm = normSpectrum(wave, data, niter=args.niter, filter_size=args.fs, trim=args.trim)

    pyfits.writeto('test.fits',
                   data_norm, header=header, clobber=True)
if __name__ == '__main__':
  main()