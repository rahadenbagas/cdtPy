# Climate Data Tools for Python v 0.2.0 (Beta Version)

A simple Python package for analyzing climate-related disciplines of Earth science data.  
Inspired by and (partially) translated from Climate Data Toolbox for MATLAB by Greene et.al., 2019 (doi:10.1029/2019GC008392).  
Disclaimer: This is a beta or pre-distribution version.

## Installation
1. Open a terminal or command prompt and navigate to the root directory of your Python package.
2. Run: pip install -e .

## Pre-requisted python libz
- numpy
- matplotlib
- xarray

## Contents
- **geoGrd (Geophysical attributes and georeferenced grids):**
  - **earth_radius** gives the nominal or latitude-dependent radius of Earth.
  - **islatlon** determines whether lat, lon is likely to represent geographical coordinates. This function is used for input parsing in many CDT functions.
  - **cdtgrid** uses meshgrid to easily create a global grid of latitudes and longitudes.
  - **cdtdim** gives the approximate dimensions of each cell in a lat,lon grid assuming a spherical Earth of radius 6371000 meters.
  - **cdtarea** gives the approximate area of each cell in a lat,lon grid assuming a spherical Earth of radius 6371000 meters. This function was designed to enable easy area-averaged weighting of large gridded climate datasets.
  - **cdtgradient** calculates the spatial gradient of gridded data equally spaced in geographic coordinates.
  - **cdtdivergence** calculates the divergence of gridded vectors on the ellipsoidal Earth's surface.
  - **cdtcurl** calculates the z component of curl for gridded vectors on an ellipsoidal Earth.
- **oceAtm (Oceans & Atmosphere):**
  - **coriolisf** returns the Coriolis frequency (also known as the Coriolis parameter or the Coriolis coefficient) for any given latitude(s).
  - **windstress** estimates wind stress on the ocean from wind speed.
  - **ekman** estimates the classical Ekman transport and upwelling/downwelling from 10 m winds.
  
(c) R.B. Hatmaja, 2023
