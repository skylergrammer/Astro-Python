#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import pyfits
import argparse
from scipy.ndimage.filters import *
from scipy.interpolate import *
from scipy.signal import correlate
import os
from scipy.constants import *

def readSpectrum(filename, wave_range, velocity, gaus_sig=2):

  # Convert km/s to m/s
  velocity = velocity * 1000  

  # Retrieve data
  data, header = pyfits.getdata(filename, header=True)
  if len(data) != int(header['naxis1']): data = data[0]

  # Smooth the input data with a 3pixel sigma guassian filter  
  data_smooth = gaussian_filter1d(data, gaus_sig)
  data_smooth_super = gaussian_filter(data, 50)

  # Convert pixels to wavelengths
  lambda1 = header['crval1']
  lambda2 = float(header['naxis1'])*header['cdelt1']+float(header['crval1'])
  
  # Convert pixels to angstroms
  raw_wave = np.arange(lambda1, lambda2, header['cdelt1'])

  # Apply a velocity
  offset = lambda v,x0: v*np.array(x0)/c  
  wave = raw_wave + offset(velocity, raw_wave)
  wave = wave[:len(data_smooth)]

  # Normalize to the flux at 5556A
  flux_5556 = data_smooth_super[np.argmin(np.abs(wave - 5556.0))]
  
  # Normalize and return only wavelengths in range
  norm_chop = [(x,y/flux_5556) for x,y in zip(wave,data_smooth) 
               if x >= wave_range[0] and x <= wave_range[1]]
  
  wave, norm_spec = zip(*norm_chop)

  return (wave, norm_spec, header)


def readStandard(filename, wave_range):  

  # Read in standard dat file
  header = ['wave', 'f1', 'sig', 'f2', 'f3']
  raw_data = np.genfromtxt(filename, usecols=(0,1,2,3,4), names=header)
  
  # Extract spectrum in wave range
  trim_data = [(x,y) for x,y,z,a,b in raw_data
               if x >= wave_range[0] and x <= wave_range[1]]  

  wave, spec = zip(*trim_data)

  # Name of file returns spectral type
  name = filename.split('/')[-1].replace('.dat','')
  
  return (wave, spec, name)


class Spectrum:
  '''
  This class keeps track of the standard spectra and their spectral type.
  '''
  def __init__(self, wave, spectrum, name):
    self.name = name
    self.wave = np.array(wave)
    self.spectrum = np.array(spectrum)


def onkey(event):
  global current
  
  # Move forwards
  if event.key == 'right' and len(all_stds) >= 2:
    fig.clf()
    if current == len(all_stds)-1: 
      current = 0
    else: 
      current += 1
    plt.plot(science.wave, science.spectrum, 'k', lw=2)
    plt.plot(all_stds[current].wave, all_stds[current].spectrum, 'r', label=all_stds[current].name)
    plt.ylim(min(science.spectrum), 3*np.mean(science.spectrum))
    plt.xlim(args.r[0], args.r[1])
    plt.legend(frameon=False, fontsize=24)    
    fig.canvas.draw()

  # Move backwards
  if event.key == 'left' and len(all_stds) >= 2:
    fig.clf()
    if current == 0: 
      current = len(all_stds)-1
    else: 
      current -= 1
    plt.plot(science.wave, science.spectrum, 'k', lw=2)
    plt.plot(all_stds[current].wave, all_stds[current].spectrum, 'r', label=all_stds[current].name)
    plt.ylim(min(science.spectrum), 3*np.mean(science.spectrum))
    plt.xlim(args.r[0], args.r[1])
    plt.legend(frameon=False, fontsize=24)
    fig.canvas.draw()
 
  # Allows user to delete models with key 'd'
  if event.key == 'd' and len(all_stds) >= 2:
    fig.clf()
    all_stds.pop(current)
    if current == len(all_stds):
      current = 0
    plt.plot(science.wave, science.spectrum, 'k', lw=2)
    plt.plot(all_stds[current].wave, all_stds[current].spectrum, 'r', label=all_stds[current].name)
    plt.ylim(min(science.spectrum), 3*np.mean(science.spectrum))
    plt.xlim(args.r[0], args.r[1])
    plt.legend(frameon=False, fontsize=24)
    fig.canvas.draw()

  # 'q' equals quit
  if event.key.lower() == 'q': exit()

###MAIN():###
parser = argparse.ArgumentParser()
parser.add_argument('--spec')
parser.add_argument('--stds', nargs='+')
parser.add_argument('--r', nargs=2, type=int, default=(3900, 6000), help='Wavelength range to fit.')
parser.add_argument('--s', type=int, default=2, help='Gaussian smoothing filter in units of Angstroms.')
parser.add_argument('--v', type=float, default=0.0, help='Velocity.')
args = parser.parse_args()

# Read in science spectrum
wave_star, spec_star, header_star = readSpectrum(args.spec, args.r, args.v, gaus_sig=args.s) 

science = Spectrum(wave_star, spec_star, header_star['object'])

# Read in all the standard spectra
all_stds = []

for each in args.stds:
  # Read standard spectrum
  try:
    wave_std, spec_std, name = readStandard(each, args.r)
    # Put data into Standard object
    std = Spectrum(wave_std, spec_std, name)
    # Append standard object to the dictionary of standards
    all_stds.append(std)
  except:
    continue

current = 0
fig,ax = plt.subplots(figsize=(15,10), dpi=72) 
fig.subplots_adjust(wspace=0.25, left=0.05, right=0.95,
                    bottom=0.125, top=0.9)
plt.xlim(args.r[0], args.r[1])
plt.ylim(min(science.spectrum), 3*np.mean(science.spectrum))
cid = fig.canvas.mpl_connect('key_press_event', onkey)
plt.plot(science.wave, science.spectrum, 'k', lw=2)
plt.plot(all_stds[current].wave, all_stds[current].spectrum, 'r', label=all_stds[current].name)
plt.legend(frameon=False, fontsize=24)
plt.show()
