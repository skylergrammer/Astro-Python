#! /usr/bin/env python

import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pyfits
import sys

def get_image(image):

  try: hdu, header = pyfits.getdata(image, header=True)
  except: sys.exit('Cannot open image.')
  
  return (hdu, header)


class Pixel:
  def __init__(self, hdu, header, center=None):
    self.hdu = hdu
    self.header = header
    self.x = self.coords(hdu, 'x')
    self.y = self.coords(hdu, 'y')
    self.z = hdu.flat[:]
    if center: self.center = center
    else: 
      self.center = (header['crpix1'], header['crpix2'])
    self.dist = self.distance()

  def coords(self, hdu, axis):
    xval, yval = np.indices(hdu.shape)
    if axis == 'x': return xval.flat[:]
    if axis == 'y': return yval.flat[:]

  def distance(self):
    x0,y0 = self.center
    r = lambda x,y,x0,y0: np.sqrt((x-x0)**2 + (y-y0)**2)
    return r(self.x, self.y, x0, y0)


def bin_dat_shit(radii, intensity, bounds=[0,20], plate_scale=1, bin_size=1):
  bin_size = bin_size*plate_scale/60.
  radii = radii*plate_scale/60.
  inner_edges = np.arange(bounds[0], bounds[1], bin_size)
  
  sky = np.mean(intensity)
  radial_curve = {'r':[], 'i':[], 'e':[]}

  for each in inner_edges:
    
    # find all the indices with radii inside r and r+bin
    vals = np.where((radii >= each) & (radii <= (each+bin_size)) & (intensity >= 0.0))
  
    mean_intensity = np.median(intensity[vals]) 
    err = np.std(intensity[vals])/np.sqrt(len(vals[0]))

    radial_curve['r'].append(each)
    radial_curve['i'].append(mean_intensity) 
    radial_curve['e'].append(err)

    # stop doing shit once SB is at sky level

  radial_curve['r'] = np.array(radial_curve['r'])
  radial_curve['i'] = np.array(radial_curve['i'])
  radial_curve['e'] = np.array(radial_curve['e'])

  return radial_curve


def fit_exponential(radial_curve):

  expfunc = lambda x,a,b,c: a*np.exp(-x/b)+c

  a0 = np.max(radial_curve['i'])
  b0 = np.max(radial_curve['r'])/5.
  c0 = 0.

 
  fit_param, cov = curve_fit(expfunc, radial_curve['r'], radial_curve['i'], 
                             p0=[a0,b0,c0])#, sigma=radial_curve['e'])
  
  fit_i = expfunc(radial_curve['r'], *fit_param)
  residuals = radial_curve['i'] - fit_i

  return (fit_i, fit_param, cov, residuals)


def main():
  
  parser = argparse.ArgumentParser()
  parser.add_argument('image', metavar='Image')
  parser.add_argument('name', type=str, default='')
  parser.add_argument('--bin', type=int, default=100., help='Bin size.  Default is 100 pixels.')
  parser.add_argument('--ps', default=1, type=float, help='Plate scale. (arcsec per pixel).')  
  parser.add_argument('--b', nargs=2, metavar='Inner Radius', type=float)
  parser.add_argument('--c', nargs=2, type=float, default=[], help='Center of radial profile, in pixels.')
  args = parser.parse_args()
  
  # Read in image
  hdu, header = get_image(args.image)

  # Put data into Pixel Class
  pix = Pixel(hdu, header, center=args.c)
  
  # Compute the radial profile
  radial_curve = bin_dat_shit(pix.dist, pix.z, bounds=args.b,
                              plate_scale=args.ps, bin_size=args.bin)

  # Fit radial profile with exponential
  fit_i, fit_param, cov, residuals = fit_exponential(radial_curve)
  print fit_param

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~PLOTS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  fig, ax = plt.subplots(figsize=(12,8), dpi=72)
  plt.subplots_adjust(wspace=0.2, left=0.1, right=0.95,
                        bottom=0.1, top=0.925)  
  
  plt.errorbar(radial_curve['r'], radial_curve['i'], yerr=radial_curve['e'],elinewidth=2,
               lw=0, marker='o', ms=5, color='k', label='Data')
  plt.plot(radial_curve['r'], fit_i, 'r', lw=3, label='Exponential Fit')

  plt.xlim(args.b)
  plt.ylim(0, np.max(radial_curve['i']))
  plt.ylabel('Median Intensity per Pixel [x$10^{-3}$]', fontsize=28)
  plt.xlabel('$\\rho$ [arcmin]', fontsize=28)
  plt.title(args.name, fontsize=28)
  plt.tick_params(axis='both', which='both',labelsize=28)
  ax.set_yticklabels(['0','5','10','15','20','25'])
  plt.show()

main()
  