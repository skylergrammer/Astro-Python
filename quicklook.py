#! /usr/bin/env python
#mod date: 2014-05-21
#upload date:
broken_mods = []

import argparse
try:
  import pyfits
except: 
  broken_mods.append('pyfits')
try: 
  import numpy as np
except: 
  broken_mods.append('numpy')
try: 
  import matplotlib.pyplot as plt
except: 
  broken_mods.append('matplotlib')
try: 
  from scipy.ndimage.filters import gaussian_filter1d
except: 
  broken_mods.append('scipy')

if broken_mods:
  print 'The following modules need to be installed: %s' % ', '.join([x for x in broken_mods])
  exit(1)


def plot_spectrum(spectrum, wave, header, color='r'):
 
  plt.plot(wave,spectrum, ls='-', color=color)
  plt.xlim(min(wave),max(wave))
  plt.ylim(min(spectrum), max(spectrum))
  plt.xlabel('$\lambda$ ($\AA$)', fontsize=24)
  #plt.ylabel('Counts', fontsize=24)
  plt.title(header['object'], fontsize=30)
 
  #plt.show()


def get_spectrum(filename):
  '''Read spectrum using pyfits.  Returns spectrum and header.
  '''
  spectrum, header = pyfits.getdata(filename, 0, header=True)

  if len(spectrum) != int(header['naxis1']): spectrum = spectrum[0]

  
  lambda1 = header['crval1']
  lambda2 = float(header['naxis1'])*header['cdelt1']+float(header['crval1'])

  wave = np.arange(lambda1, lambda2, header['cdelt1'])
  wave = wave[:len(spectrum)]
 
  return (header, wave, spectrum)


def onkey(event):
  """Allows the user to scroll through a list of spectra using the left and
     right keys.  If the current index is equal to the end of the list, it 
     simply sets the index back to 0 to allow for endless scrolling fun.  Press
     q to quit out.  It's just like IRAF's splot, just not horrible.
  """

  global current

  if event.key == 'right':
    fig.clf()
    if current == len(fileslist[:-1]): 
      current = 0
    else: 
      current += 1
    header, wave, spectrum = get_spectrum(fileslist[current])
    if args.s: spectrum = smooth(spectrum, kernel=args.s)
    plot_spectrum(spectrum, wave, header)
    fig.canvas.draw()


  if event.key == 'left':
    fig.clf()
    if current == 0: 
      current = len(fileslist[:-1])
    else: 
      current -= 1
    header, wave, spectrum = get_spectrum(fileslist[current])
    if args.s: spectrum = smooth(spectrum, kernel=args.s)
    plot_spectrum(spectrum, wave, header)
    fig.canvas.draw()
 
  if event.key == 'q': exit()


def smooth(spectrum, kernel=0):
  """Smooths the input spectrum using a user-specified Gaussian kernel.
  """

  kernel = kernel[0]
  spectrum = gaussian_filter1d(spectrum, sigma=kernel)

  return spectrum


parser = argparse.ArgumentParser(description='An easy way to view a 1D, single extension, spectrum.  Will not work on .ms.fits files.')
parser.add_argument('input', metavar='input', nargs='+', type=str, help='Filename of single spectrum or list of spectra.')
parser.add_argument('--s', nargs=1, type=str, default='3', help='Kernel size (in pix) to use for Gaussian smoothing.')
args = parser.parse_args()
fileslist = args.input

fig, ax = plt.subplots(figsize=(15,7), dpi=72) 
fig.subplots_adjust(wspace=0.25, left=0.05, right=0.95,
                    bottom=0.125, top=0.9)

current = 0
cid = fig.canvas.mpl_connect('key_press_event', onkey)  #allows key presses to be read within the plotting window.


header0, wave0, spectrum0 = get_spectrum(fileslist[0])
if args.s: spectrum0 = smooth(spectrum0, kernel=args.s)
plot_spectrum(spectrum0, wave0, header0)
plt.show()
