#! /usr/bin/env python

broke_ass_imports = []
import argparse
import sys
try:
  import numpy as np
except:
  broke_ass_imports.append('numpy')
try:
  import matplotlib.pyplot as plt
except:
  broke_ass_imports.append('matplotlib')

try:
  import pyfits
except:
  broke_ass_imports.append('pyfits')
try:
  from scipy.interpolate import *
  from scipy.ndimage.filters import *
  from scipy.stats import *
except: 
  broke_ass_imports.append('scipy')

if broke_ass_imports:
  sys.exit('\nThe following modules need to be installed to proceed: %s.\n' % (', '.join([x for x in broke_ass_imports])))
  

def getSpectrum(filename):

  # Retrieve data
  data, header = pyfits.getdata(filename, header=True)

  # Convert pixels to wavelengths
  lambda1 = header['crval1']
  lambda2 = float(header['naxis1'])*header['cdelt1']+float(header['crval1'])
  
  wave = np.arange(lambda1, lambda2, header['cdelt1'])
  wave = wave[:len(data[0])]
  
  return (wave, data[0], header)


def normSpectrum(xcoord, ycoord, header, niter=3, filter_size=50, trim=(100,100)):
  
  print '\nFilter size = %s Angstroms; N iterations = %s, and trim = %s.\n' % (filter_size, niter, trim)  

  init_size_pix = header['naxis1']

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
    
    # Gaussian filter the sobel filtered spectrum
    dyi = yi 
    dyi_smooth = gaussian_filter1d(yi, filter_size)
    
    # Fit smoothed sobel filtered spectrum with a spline
    f0 = interp1d(xi, dyi_smooth)

    # Compute the fractional difference between each data point and the spline
    frac_diff = (dyi - f0(xi))/f0(xi)

    # Compute 3sigma trimmed median and std 
    mlim = np.std(frac_diff)
    frac_diff_med = tmean(frac_diff, limits=(-3*mlim, 3*mlim))
    frac_diff_std = tstd(frac_diff, limits=(-3*mlim, 3*mlim))
    
    # Remove points that deviate from the spline by 1sigma
    xi_clean = [x for x,y,z, in zip(xi,yi,frac_diff) if np.abs(z) < frac_diff_std]
    # Regrid the data without the bad points
    yi_clean = griddata(xcoord, ycoord_sobel, xi_clean)

    # Set initial variables to the cleaned variables
    xi = xi_clean
    yi = yi_clean

    '''
    diagnostics = [(x,z) for x,y,z, in zip(xi,yi,frac_diff) if np.abs(z) > frac_diff_std]
    rej,ect = zip(*diagnostics)
    plt.plot(xi, frac_diff, 'k.')
    plt.plot([min(xi),max(xi)],[frac_diff_med,frac_diff_med],'g')
    plt.plot([min(xi),max(xi)],[-frac_diff_std,-frac_diff_std],'r')
    plt.plot([min(xi),max(xi)],[frac_diff_std,frac_diff_std],'r')
    plt.ylim(max(frac_diff), min(frac_diff))  
    plt.plot(rej, ect, 'r.')    
    plt.show()
    '''

  # Create gridded spectrum from the smoothed spectrum (non-sobel filtered)
  yi = griddata(xcoord, ycoord_smooth, xi_clean)
 
  # Crop the spectrum down to the bounds of cleaned and cropped spectrum
  ind_crop = np.where((xcoord <= max(xi)) & (xcoord >= min(xi)))  
  xcoord_crop = xcoord[ind_crop]
  ycoord_crop = ycoord[ind_crop]
 
  # Do spline interpolation of the gridded, cleaned spectrum
  f1= interp1d(xi, yi, copy=True, bounds_error=True, fill_value=np.nan)
  yn = f1(xcoord_crop)  
  
  #PLAWTS
  plt.plot(xcoord, ycoord, 'k', label='Data')
  plt.plot(xcoord, ycoord_smooth, 'g', lw=2, label='Smoothed Spectrum')
  plt.plot(xcoord_crop, yn, 'r', lw=2, label='Best Fit Continuum')
  plt.xlim(min(xcoord), max(xcoord))
  plt.ylim(0, 2*max(yn))
  plt.xlabel('$\lambda$', fontsize=30)
  plt.legend(frameon=False, fontsize=12)
  plt.title('Object:%s; N:%s; filter size:%s$\AA$' % (header['object'], niter, filter_size))
  plt.savefig('spec_'+header['object']+'_ContFit.pdf', 
              dpi=None, facecolor='w', edgecolor='w',
              transparent=True, bbox_inches=None, pad_inches=0.1,
              frameon=None)
  plt.clf()


  #Create the normalized spectrum and gaussian filter it a lil bit
  norm_spec = gaussian_filter1d(ycoord_crop/yn, 1)

  # Update the header info with correct number of pixels (some cropped off)
  header.update(NAXIS1=str(len(norm_spec)))
  # Update the header info with the correct starting wavelength
  header.update(CRVAL1=float(np.round(xcoord_crop[0], 2)))
  # Update the header with a history comment with info regarding normalization
  header.add_history('Used mk_morm_spec.py to performed normalization.')

  print '\nSpectrum went from %s pixels to %s pixels during the normalization process.\n' % (init_size_pix, header['naxis1'])
  
  return (xcoord_crop, norm_spec, header)
  
def main():

  parser = argparse.ArgumentParser(description='Take a 1D, single-extension spectrum and normalize it.')
  parser.add_argument('files', nargs='+', help='File(s) to normalize.')
  parser.add_argument('--niter', metavar='niter', type=int, default=3, help='Number of iterations to perform.')
  parser.add_argument('--fs', metavar='Size', type=int, default=100, help='Filter size in Angstroms.')
  parser.add_argument('--trim', metavar=('Left','Right'), nargs=2, type=int, default=(100,100), help='Amount to be trimmed off from the edges') 
  parser.add_argument('--clobber', action='store_true', default=False, help='If set, will clobber if necessary.')
  args = parser.parse_args()

  for each in args.files:
    # Get the un-normalized spectrum
    wave, data, header = getSpectrum(each)

    print '\n\n##########Normalizing %s##########' % (each)
    # Normalize that little bitch
    wave_norm, data_norm, header = normSpectrum(wave, data, header, niter=args.niter, filter_size=args.fs, trim=args.trim)

    # Write out normalized spectrum
    pyfits.writeto('spec_'+header['object']+'_norm.fits',
                   data_norm, header=header, clobber=args.clobber)


if __name__ == '__main__':
  main()