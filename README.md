# Ellipsoid

[![Version](https://img.shields.io/badge/Version-1.2.0-brightgreen.svg?style=plastic)](https://github.com/Robot2100/Ellipsoid/releases/tag/1.2.0)


  Calculate thermal ellipsoids and create shelx file from XDATCAR (VASP output).

#### This program uses [![Includes](https://img.shields.io/badge/Includes-1.1.6a-orange.svg)](https://github.com/Robot2100/Includes/releases/tag/1.1.6a)

### Usage:
    Ellipsoid <Filename> [-c/--cut <N>] [-q/--quiet] [-h/--help]
    
### Parameters:

  #### \<Filename\> (optional)
  File is Shelx-type file (.res/.ins). The program takes atomic coordinates and symmetry from this file.
  
  #### -c/--cut <N>
  Ignore the first steps of the "N" in XDATCAR. the default is N = 2000.
  
  #### -q/-quiet
  Print only error messages.
  
  #### -h/--help
  Print help information.
