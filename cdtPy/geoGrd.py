import numpy as np
import sys

def earth_radius(lat = None, km = False):
    """
    earth_radius() gives the radius of the Earth.

    Args:
        lat (int,float) : latitude.
        km (bool)       : units in kilometers (default units in meters)

    Synx & Desc:
        r = earth_radius() returns 6371000, the nominal radius of the Earth in meters.
        r = earth_radius(lat) gives the radius of the Earth as a function of latitude.
        r = earth_radius(..., km = True) returns values in kilometers.
    """
    # Set defaults:
    r = 6371000

    if np.any(lat):
        a = 6378137  # equatorial radius in meters
        b = 6356752  # polar radius in meters
        lat = np.radians(lat)
        r = np.sqrt(((a**2*np.cos(lat))**2 + (b**2*np.sin(lat))**2) /
            ((a*np.cos(lat))**2 + (b*np.sin(lat))**2))

    if km:
        r /= 1000

    return r

def islatlon(lat, lon):
    """
    islatlon() determines whether lat,lon is likely to represent geographical coordinates. 
    This function is use for input parsing in many CDT functions.

    Args:
        lat (int, float): latitude.
        lon (int, float): longitude.

    Synx & Desc:
        tf = islatlon(lat, lon) returns true for all input coordinates where lat
        is in the range -90 to 90 and lon is in the range -180 to 360.
    """
    tf = all([
          all(np.abs(lat[np.isfinite(lat)]) <= 90),
          all(lon[np.isfinite(lon)] <= 360),
          all(lon[np.isfinite(lon)] >= -180)])
    
    return tf

def cdtgrid(res = [1, 1], centerLon = 0):
    """
    cdtgrid() uses meshgrid to easily create a global grid of latitudes and longitudes..

    Args:
        res (int, float)        : grid resolution.
        centerLon (int, float)  : longitude as the center of the grids.

    Synx & Desc:
        [lat,lon] = cdtgrid() generates a meshgrid-style global grid of latitude and longitude
        values at a resolution of 1 degree. Postings are centered in the middle of grid cells, 
        so a 1 degree resolution grid will have latitude values of 89.5, 88.5, 87.5, etc.  
         
        [lat,lon] = cdtgrid(res) specifies grid resolution res, where res is a scalar and 
        specifies degrees. Default res is 1.
        
        [lat,lon] = cdtgrid([latres lonres]) if res is a two-element array, the first element 
        specifies latitude resolution and the second element specifies longitude resolution. 
          
        [lat,lon] = cdtgrid(...,centerLon) centers the grid on longitude value centerLon. Default
        centerLon is the Prime Meridian (0 degrees).
    """
    if np.any(res):
        assert np.max(res) < 90, 'Error, I assume. Resolution should not exceed 90 degrees.'
        if isinstance(res,(float,int)):
            res = [res, res]
        else:
            res = res

    # Generate grid:
    lat1 = np.arange(90 - res[0] / 2, -90 - res[0] / 2, -res[0])
    lon1 = np.arange(-180 + res[1] / 2, 180 - res[1] / 2 + res[1], res[1])
    lat, lon = np.meshgrid(lat1, lon1)
    
    # Recenter:
    if centerLon != 0:
        lon += centerLon
        
    return lat, lon

def cdtdim(lat, lon, km = False):
    """
    cdtdim() gives the approximate dimensions of each cell in a lat,lon grid assuming an
    ellipsoidal earth.
    
    Args:
        lat (int, float): latitude.
        lon (int, float): longitude.
        km (bool)       : units in kilometers (default units in meters)

    Synx & Desc:
        [dx, dy] = cdtdim(lat, lon) gives an approximate dimensions in meters of each grid cell given 
        by the geographical coordinate grids lat, lon. Inputs lat and lon must have matching dimensions,
        as if they were created by meshgrid. 
        [dx, dy] = cdtdim(lat, lon, km = True) gives grid cell sizes in kilometers rather than the default meters.
    """
    # Initial errors check:
    assert len(sys.argv) >= 2, 'Error: cdtgradient requires at least two inputs: latNo documeNo docume and lon.'
    assert np.array_equal(np.ndim(lat),2) and np.array_equal(np.ndim(lon),2), 'Error: lat and lon must be grids, not arrays.'
    assert np.array_equal(np.shape(lat),np.shape(lon)), 'Input error: the dimensions of lat and lon must match.'
    assert islatlon(lat, lon) == True, 'Input error: cdtgradient requires the first two inputs to be lat and lon.'
    
    # Set defaults
    r = earth_radius(lat)
    
    # Input parsing
    if km:
        r /= 1000
        
    # Determine grid sizes dlat and dlon: 
    
    [dlatY, dlatX] = np.gradient(lat)
    [dlonY, dlonX] = np.gradient(lon) 
    
    # We don't know if lat and lon were created by [lat,lon] = meshgrid(lat_array,lon_array) or [lon,lat] = meshgrid(lon_array,lat_array) 
    # but we can find out: 
    if np.array_equal(dlatX, np.zeros_like(lat)):
        dlat = dlatY
        dlon = dlonX
        assert np.array_equal(dlonY, np.zeros_like(lon)), 'Error: lat and lon must be monotonic grids, as if created by np.meshgrid.'
    else:
        dlat = dlatX
        dlon = dlonY
        assert np.array_equal(dlonX, np.zeros_like(lon)) or np.array_equal(dlatY, np.zeros_like(lon)), 'Error: lat and lon must be monotonic grids, as if created by np.meshgrid.'
    
    # Calculate dimensions based on dlat and dlon: 
    
    dy = dlat * r * np.pi / 180
    dx = (dlon / 180) * np.pi * r * np.cos(np.deg2rad(lat))
    
    return dx, dy

def cdtarea(lat, lon, km2 = False):
    """
    cdtarea() gives the approximate area of each cell in a lat,lon grid assuming an ellipsoidal
    earth. This function was designed to enable easy area-averaged weighting of large gridded 
    climate datasets.
    
    Args:
        lat (int, float): latitude.
        lon (int, float): longitude.
        km2 (bool)      : units in square kilometers (default units in square meters)

    Synx & Desc:
        A = cdtarea(lat, lon) gives an approximate area of each grid cell given by lat,lon. Inputs
        lat and lon must have matching dimensions, as if they were created by meshgrid. 
        A = cdtarea(lat, lon , km2 = True) gives grid cell area in square kilometers.
    """
    
    [dx, dy] = cdtdim(lat, lon, km = km2)
    
    A = np.abs(dx * dy)

    return A

def cdtgradient(lat, lon, F, km = False):
    """
    cdtgradient calculates the spatial gradient of gridded data equally spaced in geographic
    coordinates. 
    
    Args:
        lat (int, float): latitude.
        lon (int, float): longitude.
        F (int, float)  : gridded variable F.
        km (bool)       : units in kilometers (default units in meters)

    Synx & Desc:
        [FX, FY] = cdtgradient(lat, lon, F) for the gridded variable F and corresponding geographic
        coordinates lat and lon, cdtgradient calculates FX, the spatial rate of west-east change 
        in F per meter along the Earth's surface, and FY, the south-north change in F  per meter. 
        This function assumes an ellipsoidal Earth as modeled by the cdtdim and earth_radius. A
        positive value of FX indicates that F increases from west to east in that grid cell, as
        a positive value of FY indicates F increases from south to north. F can be a 2D or 3D matrix
        whose first two dimensions must correspond to lat and lon. If F is 3D, outputs FX and FY
        will also be 3D, with each grid along the third dimension calculated separately. 

        [FX, FY] = cdtgradient(lat, lon, F, km = True) returns gradients per kilometer rather than the
        default meters. 
    """
    # Initial errors check: 
    assert len(sys.argv) >= 2, 'Error: cdtgradient requires at least three inputs: lat, lon, and F.'
    assert len(sys.argv) <= 4, 'Unrecognized input.'
    assert np.array_equal(np.ndim(lat),2) and np.array_equal(np.ndim(lon),2), 'Error: lat and lon must be grids, not arrays.'
    assert islatlon(lat, lon) == True, 'Input error: cdtgradient requires the first two inputs to be lat and lon.'
    
    # Estimate grid size: 
    if km:
        [dx, dy] = cdtdim(lat, lon, km = True)
    else:
        [dx, dy] = cdtdim(lat, lon)

    # Calculate gradient: 
    if F.ndim == 3:
        assert np.array_equal(np.shape(lat),np.shape(F[0,:,:])) or np.array_equal(np.shape(lon),np.shape(F[0,:,:])), 'Input error: dimensions of lat, lon, and F must match.'
        [_,FY,FX] = np.gradient(F)
    else:    
        [FX, FY] = np.gradient(F)
        
    FX = FX/dx
    FY = FY/dy

    return FX, FY

def cdtdivergence(lat, lon, U, V):
    """
    cdtdivergence() calculates the divergence of gridded vectors on the ellipsoidal
    Earth's surface. 
    
    Args:
        lat (int, float): latitude.
        lon (int, float): longitude.
        U (int, float)  : gridded vector's component U.
        V (int, float)  : gridded vector's component V.

    Synx & Desc:
        D = cdtdivergence(lat,lon,U,V) uses cdtdim to estimate the dimensions of each
        grid cell in the |lat,lon| grid, then computes the divergence of the gridded 
        vectors |U,V|. Units of D are the units of U and V divided by meters. 
    """
    # Initial error checks: 
    assert len(sys.argv) >= 4, 'Error: cdtdivergence only requires four inputs: lat, lon, U, and V.'
    assert np.array_equal(np.ndim(lat),2) and np.array_equal(np.ndim(lon),2), 'Error: lat and lon must be grids, not arrays.'
    assert np.array_equal(np.shape(lat),np.shape(lon)), 'Input error: the dimensions of lat and lon must match.'
    assert np.array_equal(np.shape(U),np.shape(V)), 'Input error: the dimensions of U and V must match.'
    assert islatlon(lat, lon) == True, 'Input error: Some of the values in lat or lon do not match typical lat,lon ranges. Check inputs and try again.'

    # Estimate grid size:
    [dx, dy] = cdtdim(lat, lon)
    
    # Calculate:
    if U.ndim == 3 and V.ndim == 3:
        assert np.array_equal(np.shape(lat),np.shape(U[0,:,:])) or np.array_equal(np.shape(lon),np.shape(U[0,:,:])), 'Input error: dimensions of lat, lon, U, and V must match.'
        [_,_,UX] = np.gradient(U)
        [_,VY,_] = np.gradient(V)
    else:    
        [UX,_] = np.gradient(U)
        [_,VY] = np.gradient(V)
        
    UX = UX/dx
    VY = VY/dy
    D = UX + VY
    
    return D

def cdtcurl(lat, lon, U, V):
    """
    cdtcurl() calculates the z component of curl for gridded vectors on an ellipsoidal Earth. 
    
    Args:
        lat (int, float): latitude.
        lon (int, float): longitude.
        U (int, float)  : gridded vector's component U.
        V (int, float)  : gridded vector's component V.

    Synx & Desc:
        Cz = cdtcurl(lat,lon,U,V) uses cdtdim to estimate the dimensions of each grid cell in the 
        lat,lon grid, then computes the curl of the gridded vectors U,V.  
    """
    # Initial errors check: 
    assert len(sys.argv) >= 4, 'Error: cdtdivergence only requires four inputs: lat, lon, U, and V.'
    assert np.array_equal(np.ndim(lat),2) and np.array_equal(np.ndim(lon),2), 'Error: lat and lon must be grids, not arrays.'
    assert np.array_equal(np.shape(lat),np.shape(lon)), 'Input error: the dimensions of lat and lon must match.'
    assert np.array_equal(np.shape(U),np.shape(V)), 'Input error: the dimensions of U and V must match.'
    assert islatlon(lat, lon) == True, 'Input error: Some of the values in lat or lon do not match typical lat,lon ranges. Check inputs and try again.'

    # Estimate grid size:
    [dx, dy] = cdtdim(lat, lon)
    
    # Calculate:
    if U.ndim == 3 and V.ndim == 3:
        assert np.array_equal(np.shape(lat),np.shape(U[0,:,:])) or np.array_equal(np.shape(lon),np.shape(U[0,:,:])), 'Input error: dimensions of lat, lon, U, and V must match.'
        [_,UY,_] = np.gradient(U)
        [_,_,VX] = np.gradient(V)
    else:    
        [_,UY] = np.gradient(U)
        [VX,_] = np.gradient(V)
        
    UY = UY/dy
    VX = VX/dx
    Cz = VX - UY
    
    return Cz