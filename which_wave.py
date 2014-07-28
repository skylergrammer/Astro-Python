#! /usr/bin/env python

import numpy as np
import argparse


def mkRestWave(wavelength, velocity):

  rest_wave = lambda x,v: x / (1. + v/2.98e5)

  return rest_wave(wavelength, velocity)


def mkObsWave(wavelength, velocity):

  obs_wave = lambda x,v:  (1. + v/2.98e5)*x

  return obs_wave(wavelength, velocity)


def main():

  parser = argparse.ArgumentParser(description='Returns wavelengths that have been blue or red shifted by a given wavelength.')
  parser.add_argument('wave', nargs='+', type=float, help='Wavelength(s), in Angstroms, to be converted to rest or observed given a velocity in km/s.')
  parser.add_argument('--v', type=float, default=0, help='Velocity (km/s).')
  parser.add_argument('--o', action='store_true', help='Set if using observed wavelengths.')
  parser.add_argument('--r', action='store_true', help='Set if using rest wavelengths.')

  args = parser.parse_args()

  output = []
  for each in args.wave:
    
    if args.o:
      rest_wave = mkRestWave(each, args.v)
      output.append(rest_wave)

    if args.r:
      obs_wave = mkObsWave(each, args.v)
      output.append(obs_wave)

  print '\nShifted wavelength(s) given a velocity of %s km/s: \n%s\n' % (args.v, ', '.join(["%.2f" % x for x in output]))

if __name__ == '__main__':
  main()
