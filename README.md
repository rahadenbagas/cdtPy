# Climate Data Tools for Python (cdtPy) v 1.0.0

A simple Python package for analyzing climate-related disciplines of Earth science data. 
Inspired by and partially translated from Climate Data Toolbox for MATLAB by Greene et.al., 2019 (doi:10.1029/2019GC008392). 
Unravel the mysteries of our planet's climate with cdtPy. 
Harness the power of data analysis, visualization, and research to contribute to our collective understanding of Earth's intricate climate systems. 

Why Choose **cdtPy**:

**User-Friendly**: Our user-centric design ensures that both beginners and seasoned professionals can navigate the package effortlessly.
**Open Source**: cdtPy is an open-source project, promoting collaboration and constant improvement within the Earth science community.
**Robust Documentation**: Our extensive documentation, tutorials, and examples guide you through every aspect of using the package, from installation to advanced analysis.

## Installation
1. pip install cdtPy

## Pre-requisted python libz
- numpy
- matplotlib
- xarray

## Unleash the full potential of cdtPy:
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
