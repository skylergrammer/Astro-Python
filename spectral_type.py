#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import pyfits
import argparse
from scipy.ndimage.filters import *
from scipy.interpolate import *
from PyAstronomy import pyasl
from scipy.signal import correlate
import os


def readSpectrum(filename, wave_range):

  # Retrieve data
  data, header = pyfits.getdata(filename, header=True)
  
  # Smooth the input data with a 3pixel sigma guassian filter  
  data_smooth = gaussian_filter1d(data[0], 5)

  # Convert pixels to wavelengths
  lambda1 = header['crval1']
  lambda2 = float(header['naxis1'])*header['cdelt1']+float(header['crval1'])
  
  # Convert pixels to angstroms
  wave = np.arange(lambda1, lambda2, header['cdelt1'])
  wave = wave[:len(data_smooth)]

  # Normalize to the flux at 5556A
  flux_5556 = data_smooth[np.argmin(np.abs(wave - 5556.0))]
  
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
  name = filename.split('/')[-1].strip('.dat')
  
  return (wave, spec, name)


class Spectrum:
  '''
  This class keeps track of the standard spectra and their spectral type.
  '''
  def __init__(self, wave, spectrum, name):
    self.name = name
    self.wave = np.array(wave)
    self.spectrum = np.array(spectrum)


def xcorrelate(science, standard):

  new_wave = np.linspace(min(science.wave), max(science.wave), num=len(standard.wave))

  new_spec = griddata(science.wave, science.spectrum, list(new_wave))


def main():

  parser = argparse.ArgumentParser()
  
  parser.add_argument('spectrum')
  parser.add_argument('standards', nargs='+')
  parser.add_argument('--r', nargs=2, type=int, required=True)
  
  args = parser.parse_args()
  
  # Read in science spectrum
  wave_star, spec_star, header_star = readSpectrum(args.spectrum, args.r) 
  
  science = Spectrum(wave_star, spec_star, header_star['object'])


  # Read in all the standard spectra
  all_stds = {}
  
  for each in args.standards:
    # Read standard spectrum
    wave_std, spec_std, name = readStandard(each, args.r)

    # Put data into Standard object
    std = Spectrum(wave_std, spec_std, name)
   
    # Append standard object to the dictionary of standards
    all_stds[std.name] = std

  
  for each in all_stds: 
    fig,ax = plt.subplots(figsize=(15,10), dpi=72) 
    fig.subplots_adjust(wspace=0.25, left=0.05, right=0.95,
                        bottom=0.125, top=0.9)
    plt.plot(science.wave, science.spectrum, 'k', lw=2)
    plt.plot(all_stds[each].wave, all_stds[each].spectrum, 'r', label=each)
    plt.legend()
    plt.show()
    output = xcorrelate(science, all_stds[each])


if __name__ == '__main__':
  main()