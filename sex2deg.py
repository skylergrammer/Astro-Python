#! /usr/bin/env python

def sex2deg(coord, RA=True):

  if RA: 
    rah, ram, ras = coord.split(':')
    deg = (float(rah) + float(ram)/60. + float(ras)/3600.) * 15.
  if not RA:
    ded, dem, des = coord.split(':')
    deg = float(ded) + float(dem)/60. + float(des)/3600.
  
  return float(deg)
