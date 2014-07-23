#! /usr/bin/env python
import argparse
import pyfits
import numpy as np
import os
import matplotlib.pyplot as plt
from scipy.ndimage.filters import gaussian_filter1d


def readSpectrum(filename, smooth=False):
    '''
    Read in a spectrum.  Return wavelengths, intensities, and header.
    '''
    data, header = pyfits.getdata(filename, 0, header=True)
    if len(data) != int(header['naxis1']):
        data = data[0]

    waves = np.linspace(header['crval1'],header['naxis1']*header['cdelt1']+header['CRVAL1'], num=header['naxis1'])

    if smooth:
      data = gaussian_filter1d(data, sigma=3.0)
      print 'Smoothed'

    return (waves, data, header)


def main():
    parser = argparse.ArgumentParser(description='Copy first extension of input files to target directory.')
    
    parser.add_argument('filelist',nargs='+',help='List of files with 1D spectra.')
    parser.add_argument('--s', action='store_true', default=False, help='If set, does a 3 pixel Gaussian smooth.')
    parser.add_argument('-outdir',default='.',help="Target directory to copy files (default = './')")


    args = parser.parse_args()

    fcounts = 0
    for filename in args.filelist:
        waves,data,hdr = readSpectrum(filename, smooth=args.s)
        
        #outname
        outname = os.path.splitext(filename)[0] + '.txt'
        outname = os.path.join(args.outdir,outname)

        print 'Writing to %s' % outname
        f = open(outname,'w')
        for x,y in zip(waves,data):
            f.write('%f  %f\n' % (x,y))
        f.close()

        fcounts += 1

    print '%i spectra written' % fcounts

if __name__ == '__main__':
    main()
