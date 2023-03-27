# Climate Data Tools for Python v 0.1.0

A simple Python package for analyzing of climate-related disciplines of Earth science data.

## Installation
pip install cdtPy

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

  
(c) R.B. Hatmaja, 2023
