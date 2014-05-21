#! /usr/bin/env python

import argparse
import sys
try:
  import numpy as np
except:
  sys.exit("Error: Please install numpy to proceedl.")


def get_data(filename, columns=[1,2,3]):

  try:
    data = np.genfromtxt(filename, names=['ID', 'RA', 'DEC'], dtype=['S30', 'S30', 'S30'], usecols=columns, skip_header=True)

  except:
    sys.exit("Error: There was a problem reading in your file.")

  return data
  


def print_regions(data, coords='fk5', shape='circle', size='0.25"'):

  with open("temp.reg",'w+') as f:
    
    f.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
    f.write(coords+'\n')
    
    for targid,ra,dec in data:
      meat = ','.join(['('+ra,dec,size+')'])
      end = ''.join([' # text={',targid,'}'])
    
      f.write(''.join([shape,meat,end,'\n']))

  return

def main():

  parser = argparse.ArgumentParser()
  parser.add_argument('f', metavar='Filename', help="File which contains targets.")
  parser.add_argument('--col', metavar='Column', nargs=3, type=int, default=[0,1,2], help="Columns which contain the ID, x, and y for each target. Default = [0,1,2] which corresponds to the first, second, and third columns.")
  parser.add_argument('--type', required=True, choices=('image','fk5'))
  
  args = parser.parse_args()

  data = get_data(args.f, columns=args.col)
  print_regions(data, shape='circle', size='0.25"')


if __name__ == '__main__':
  main()
