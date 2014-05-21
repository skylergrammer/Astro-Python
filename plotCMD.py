#! /usr/bin/env python
#mod-date:2014-04-17
#upload-date: 2014-04-17
#version 2.1

import sys
import matplotlib.pyplot as plt
import numpy as np
import pyfits
import argparse
import matplotlib as mpl
mpl.rcParams['text.usetex']=True
mpl.rcParams['text.latex.unicode']=True
from matplotlib.ticker import MultipleLocator, MaxNLocator

def read_cmd(filename):

  data, header = pyfits.getdata(filename, header=True)

  return (data, header)


def plot_cmd(data, header, show='screen'):

  targ = np.where(data['hstid'] == header['hstid'])
  
  fig, ax = plt.subplots(nrows=1, ncols=2, sharex=True, sharey=True, figsize=(12,12), dpi=72)
  fig.subplots_adjust(wspace=0.2, left=0.08, right=0.95,
                        bottom=0.1, top=0.925)

  ax[0].plot(data['BV'], data['V'], 'ko', ms=4, aa=True)
  ax[0].plot(data['BV'][targ], data['V'][targ], color='r', marker='o', ms=20, alpha=0.5)
  ax[0].errorbar(data['BV'], data['V'], yerr=data['Verr'], xerr=data['BVerr'], aa=True,
                 elinewidth=2, capsize=1, ls='none', color='k')  

  ax[1].plot(data['VI'], data['V'], 'ko', ms=4, aa=True)
  ax[1].plot(data['VI'][targ], data['V'][targ], color='r', marker='o', ms=20, alpha=0.5)
  ax[1].errorbar(data['VI'], data['V'], yerr=data['Verr'], xerr=data['VIerr'], aa=True,
                 elinewidth=2, capsize=1, ls='none', color='k') 

  fig.subplots_adjust(wspace=0)

  # Set axis limits
  ax[0].set_ylim(26, 18)
  ax[0].set_xlim(-1, 3)
  ax[1].set_xlim(-1, 3)  
  
  ax[0].tick_params(axis='both', which='both', labelsize=28)
  ax[1].tick_params(axis='y', which='both', labelleft='off')
  ax[1].tick_params(axis='x', which='both', labelsize=28)

  # Set the major and minor axis tick widths
  xmajorLocator = MultipleLocator(1)
  xminorLocator = MultipleLocator(0.25)
  ymajorLocator = MultipleLocator(1)
  yminorLocator = MultipleLocator(0.5)
  ax[0].xaxis.set_major_locator(MaxNLocator(integer=True, prune='lower'))
  ax[0].yaxis.set_minor_locator(yminorLocator)


  # Set the x and y labels
  ax[0].set_ylabel('$V$', fontsize=36)
  ax[0].set_xlabel('$(B-V)$', fontsize=36)
  ax[1].set_xlabel('$(V-I)$', fontsize=36)
  plt.suptitle(header['hstid'], fontsize=36)  

  # Save the figure

  if show != 'screen': plt.savefig(header['lbtid']+'_cmd.'+ show, 
                    dpi=None, facecolor='w', edgecolor='w',
                    transparent=True, bbox_inches=None, pad_inches=0.1,
                    frameon=None)

  else: plt.show()


def main():

  parser = argparse.ArgumentParser(description="Creates CMDs from fits tables of optical photometry.  plotCMD.py requires the following columns: HSTID LBTID Vmag B-V V-I.")
  parser.add_argument('cmd', nargs='+', help="Fits table containing the data.")
  parser.add_argument('--show', type=str, choices=['pdf', 'png', 'jpg', 'screen'], default='screen', help="The output file format. Default is a pdf.")
  
  args = parser.parse_args()

  for each in args.cmd:
    data, header = read_cmd(each)
    plot_cmd(data, header, show=args.show)


if __name__ == '__main__':
  main()
