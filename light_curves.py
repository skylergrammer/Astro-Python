#!/usr/bin/env python
#mod date: 2014-04-17
#upload date: 2014-04-17
#version: 1.2
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import argparse


def get_ref_data(filename):
  """Reads the file that contains the reference information for each object -- e.g. id, light curve filenames, reference mags, etc.
  """
  header = ['x_shifted','y_shifted','id',   #the header column names
            'lcR','xR','yR','magR','emagR',
            'lcV','xV','yV','magV','emagV',
            'lcB','xB','yB','magB','emagB',
            'lcU','xU','yU','magU','emagU']

  fmt = ['f8','f8','S30', #format keys for each column
         'S40','f8','f8','f8','f8',
         'S40','f8','f8','f8','f8',
         'S40','f8','f8','f8','f8',
         'S40','f8','f8','f8','f8',]

  key = np.genfromtxt(filename, names=header, autostrip=True, dtype=fmt)

  return key


def gen_lc(ref_data, sourceid, zps, dm, show="screen"):
  """Given the reference data, will find the sourceid and locate all the files containing the UBVR light curves.  Feeds the light curve to plot_lc to do the actual plotting.
  """
  #If distance modulus not set, set it to that of M101
  if not dm: dm = 29.05
  
  #Check to see if provided ID is actually in the reference file
  all_ids = [x['id'] for x in ref_data]
  if not sourceid in all_ids:
    sys.exit('%s not in reference file' %sourceid)

  #Get the reference data corresponding to provided ID
  ref_line = [x for x in ref_data if sourceid in x]
  ref_line = ref_line[0]
  
  if os.path.isfile(ref_line['lcU']):
    header = ['mjd', 'dCounts', 'edCounts']
    lc = np.genfromtxt(ref_line['lcU'], usecols=[0,1,2], names=header)
    didplotU = plot_lc((ref_line['magU'], ref_line['emagU']), sourceid, lc, dm, 
              (zps['U'],zps['eU']), bandpass='U')

    if show != "screen" and didplotU: 
      plt.savefig(sourceid+'_lcU.'+show.lower(), 
                  dpi=None, facecolor='w', edgecolor='w',
                  transparent=True, bbox_inches=None, pad_inches=0.1,
                  frameon=None)

  if os.path.isfile(ref_line['lcB']): 
    header = ['mjd', 'dCounts', 'edCounts']
    lc = np.genfromtxt(ref_line['lcB'], usecols=[0,1,2], names=header)
    didplotB = plot_lc((ref_line['magB'], ref_line['emagB']), sourceid, lc, dm,
            (zps['B'],zps['eB']), bandpass='B')

    if show != "screen" and didplotB: 
      plt.savefig(sourceid+'_lcB.'+show.lower(), 
                  dpi=None, facecolor='w', edgecolor='w',
                  transparent=True, bbox_inches=None, pad_inches=0.1,
                  frameon=None)   

  if os.path.isfile(ref_line['lcV']): 
    header = ['mjd', 'dCounts', 'edCounts']
    lc = np.genfromtxt(ref_line['lcV'], usecols=[0,1,2], names=header)
    didplotV = plot_lc((ref_line['magV'], ref_line['emagV']), sourceid, lc, dm, 
            (zps['V'],zps['eV']), bandpass='V')

    if show != "screen" and didplotV: 
      plt.savefig(sourceid+'_lcV.'+show.lower(), 
                  dpi=None, facecolor='w', edgecolor='w',
                  transparent=True, bbox_inches=None, pad_inches=0.1,
                  frameon=None)

  if os.path.isfile(ref_line['lcR']): 
    header = ['mjd', 'dCounts', 'edCounts']
    lc = np.genfromtxt(ref_line['lcR'], usecols=[0,1,2], names=header)
    didplotR = plot_lc((ref_line['magR'], ref_line['emagR']), sourceid, lc, dm,
            (zps['R'],zps['eR']), bandpass='R')

    if show != "screen" and didplotR: 
      plt.savefig(sourceid+'_lcR.'+show.lower(), 
                  dpi=None, facecolor='w', edgecolor='w',
                  transparent=True, bbox_inches=None, pad_inches=0.1,
                  frameon=None)
  
  if show == 'screen': plt.show()


def plot_lc(ref_info, sourceid, lc, dm, zp_info, bandpass='V'):
  """From the provided light curve, zero points, and distance modulus, plot the light curve.
  """

  #set plot colors depending on bandpass
  if bandpass == 'U': col='#9E00FF'
  elif bandpass == 'B': col='#0046FF'
  elif bandpass == 'V': col='#1DBD00'
  elif bandpass == 'R': col='#FF4500'
  
  ref_mag, ref_mag_err = ref_info
  zp, ezp = zp_info
  offset = zp + 25.0
  
  #convert ref mag to ref counts
  ref_counts = 10**((ref_mag-offset + dm)/-2.5)
  
  #convert dcounts to a magnitude
  counts_i = [ref_counts - x if x < ref_counts else np.nan 
              for x in lc['dCounts']]
  mag_i = np.array([-2.5*np.log10(x)+offset for x in counts_i])
  Beta = (2.5*ref_mag_err)**2
  emag_i = [np.sqrt(Beta + (1.1*y/x)**2 + ezp**2) 
            for x,y in zip(counts_i,lc['edCounts'])]
   
  rms = lambda x,x0: np.sqrt(np.sum(np.power(x-x0, 2))/len(x))
  mean_mag = np.median(mag_i)
  mean_error = np.median(emag_i)  
  rms_mag = rms(mag_i, mean_mag)
  max_dev = np.abs(max(mag_i) - min(mag_i))

  # If the RMS magnitude is greater than the average photometric error, plot
  if rms_mag >= mean_error:

    #major x-axis interval set to 365 days
    xmajorLocator = MultipleLocator(365)
    #major y-axis interval set to 0.1 mags
    #ymajorLocator = MultipleLocator(np.round(max_dev/2.5,1))
    #minor y-axis interval set to 0.05 mags
    #yminorLocator = MultipleLocator(np.round(max_dev/5,2))


    fig, ax = plt.subplots(figsize=(14,8), dpi=72)
    fig.subplots_adjust(wspace=0, left=0.15, right=0.98,
                          bottom=0.15, top=0.95)

    #plot lightcurve onto a grid
    plt.plot(lc['mjd'], mag_i, 'ko', ms=15,
             ls='--', lw=5, label=sourceid)
    plt.legend(frameon=False, fontsize=30)
    plt.errorbar(lc['mjd'], mag_i, yerr=emag_i, elinewidth=5, capsize=0, color='k', alpha=0.5)  
    plt.grid(b=True, which='major', ls='-')
    plt.grid(b=True, which='minor', ls=':')

    #set y limits to mean magnitude plus or minus 0.25 mags
    plt.ylim(np.round(max(mag_i)+3*mean_error, decimals=2), 
             np.round(min(mag_i)-3*mean_error, decimals=2))

    plt.xlabel('MJD [days]', fontsize=36)
    ax.yaxis.labelpad = 10
    plt.ylabel(' '.join([bandpass,'[mag]']),fontsize=36)

    #Set major and minor axis intervals
    ax.xaxis.set_major_locator(xmajorLocator)
    #ax.yaxis.set_major_locator(ymajorLocator)
    #ax.yaxis.set_minor_locator(yminorLocator)

    ax.tick_params(axis='both', which='both', labelsize=28)
    ax.ticklabel_format(style='plain', axis='both', useOffset=False)
    #plt.suptitle(sourceid, fontsize=30)  

    return True
  
  else:
    return False
  
def main():

  parser = argparse.ArgumentParser()
  parser.add_argument('ref', metavar='reference_file')
  parser.add_argument('--chip', choices=['1','2','3','4'], required=True, help='The chip number.  This is important for determining the photometric zeropoints.')
  parser.add_argument('--id', type=str,  help='Options: single source id, "all", or "list".  "all" will make light curves for all the sources and "list" will list all the sources, then quit".')
  parser.add_argument('--dm', type=float, help='Distance modulus.  Default is set to 29.05 (M101).')  
  parser.add_argument('--show', choices=["screen", "pdf", "png", "jpg"], default="screen" , help='How the light curve will be viewed.')

  args = parser.parse_args()

  #The zeropoints for each filter
  if args.chip == '1':  
    zps = {'U':7.02, 'eU':0.01, 
           'B':8.65, 'eB':0.01,
           'V':9.22, 'eV':0.01,
           'R':8.74, 'eR':0.01}
  if args.chip == '2':  
    zps = {'U':6.95, 'eU':0.01, 
           'B':9.43, 'eB':0.07,
           'V':9.42, 'eV':0.01,
           'R':8.94, 'eR':0.03}
  if args.chip == '3':  
    zps = {'U':7.00, 'eU':0.01, 
           'B':8.62, 'eB':0.01,
           'V':9.11, 'eV':0.01,
           'R':8.74, 'eR':0.01}
  if args.chip == '4':  
    zps = {'U':7.02, 'eU':0.01, 
           'B':8.67, 'eB':0.01,
           'V':9.17, 'eV':0.01,
           'R':8.78, 'eR':0.01}

  #Retrieves the reference data
  ref_data =  get_ref_data(args.ref)  #read the ref file

  #Lists all the source IDs in the reference data
  if args.id.lower() == "list": 
    for line in ref_data['id']: print line
    exit()

  #Generate light curves
  if args.id.lower() == 'all':
    for line in ref_data['id']: 
      gen_lc(ref_data, line, zps, args.dm, show=args.show)
  else: 
    gen_lc(ref_data, args.id, zps, args.dm, show=args.show)

if __name__ == '__main__':
  main()
