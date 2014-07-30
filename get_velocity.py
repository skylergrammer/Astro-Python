#! /usr/bin/env python

import numpy as np
import argparse


def getVel(obs_wave, rest_wave):

  v = lambda x,x0: ((x-x0)/x0) * 2.98e5

  return v(obs_wave, rest_wave)


def main():

  parser = argparse.ArgumentParser(description='Returns wavelengths that have been blue or red shifted by a given wavelength.')
  parser.add_argument('--o', metavar='Observed', type=float, required=True)
  parser.add_argument('--r', metavar='Rest', type=float, required=True)
  args = parser.parse_args()

  output = getVel(args.o, args.r)
  velocity = str(output)
  print '\nVelocity is %s km/s.\n' % velocity

if __name__ == '__main__':
  main()
