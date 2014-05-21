#! /usr/bin/env python

import argparse
import pyfits

def get_data(filename):

  data = pyfits.getdata(filename)

  return data


def mk_guide_file(data):

  f = open('hecto_guides.dat', 'w+')
  for i,line in enumerate(data):
    ra = ':'.join(line['_RAJ2000'].split())
    dec = ':'.join(line['_DEJ2000'].split())
    object = 'guide'+str(i+1).format('S5')
  
    foo = '\t\t'.join([ra, dec, object, '--', line['Jmag'].astype('|S5'), 'Guide',])
    f.write(foo+'\n')

  f.close()

def main():

  parser = argparse.ArgumentParser()
  parser.add_argument('filename')
  args = parser.parse_args()

  data = get_data(args.filename)

  mk_guide_file(data)

main()