import numpy as np
import sys
import cdtPy.geoGrd as geoGrd

def coriolisf(lat, rot = 7.2921e-5):
    """
    coriolisf() returns the Coriolis frequency for any given latitude(s). 
    The Coriolis frequency is sometimes called the Coriolis parameter or the 
    Coriolis coefficient.

    Args:
        lat (int, float): latitude.
        rot (int, float): Earth's rate of rotation.

    Synx & Desc:
        f = coriolisf(lat) returns Earth's Coriolis frequency for any locations(s) 
        at latitude(s) lat. By default, the Coriolis frequency is given in units of
        rad/s. 

        f = coriolisf(lat, rot) specifies a rate of rotation rot. By default, rot
        is Earth's present-day rate of rotation 7.2921 x 10^-5 rad/s, but a different
        rate may be specified to model Earth at a different time, other celestial 
        bodies, or the Coriolis parameter in a different unit of frequency.
    """
    # Initial error checks:
    assert np.max(np.abs(lat)) <= 90, 'Latitude value(s) out of realistic bounds. Check inputs and try again.'

    # Input parsing:
    if len(sys.argv) > 1:
        assert isinstance(rot,(float,int)), 'Error: rotation rate rot must be numeric.'

	# Coriolis calculation:
    lat = np.radians(lat)
    f = 2*rot*np.sin(lat)
    
    return f

def windstress(u10, v10, Cd = 1.25e-3, rhoa = 1.225):
    """
    windstress() estimates wind stress on the ocean from wind speed.

    Args:
        u10 (int, float): zonal 10-m wind component.
        v10 (int, float): meridional 10-m wind component.
        Cd (int, float) : drag coefficient of friction, default Cd is 1.25e-3.
        rhoa (int, float): air density, default rho is 1.225 kg/m**3

    Synx & Desc:
        [Taux,Tauy] = windstress(u10, v10) simultaneously computes zonal and meridional components of wind stress.

        [...] = windstress(...,Cd = Cd) specifies a coefficient of friction Cd. Default Cd is 1.25e-3, which is a 
        global average (Kara et al., 2007) but in reality Cd can vary quite a bit in space and time. Cd can be a 
        scalar or a vector, 2D matrix, or 3D matrix the same size as u10 (and v10 if v10 is included). 

        [...] = windstress(...,rhoa = rhoa) specifies air density, which can be a scalar or a vector, 2D matrix, 
        or 3D matrix the same size as u10 (and v10 if v10 is included). Default value of rho is 1.225 kg/m**3.
    """
    # Calculate wind stress
    U = np.hypot(u10, v10)
    Tau = rhoa * Cd * (u10**2 + v10**2)
    
    # TauX
    TauX = Tau * u10 / U
    # TauY
    TauY = Tau * v10 / U
    
    return TauX, TauY

def ekman(lat, lon, u10, v10, Tau = False, Cd = 1.25e-3, rhow = 1025):
    """
    ekman() estimates the classical Ekman transport and upwelling/downwelling from 10 m winds.

    Args:
        lat (int, float): latitude.
        lon (int, float): longitude.
        u10 (int, float): zonal 10-m wind component.
        v10 (int, float): meridional 10-m wind component.
        Tau (bool)      : if wind stress data are used, instead of 10-m wind.
        Cd (int, float) : drag coefficient of friction, default Cd is 1.25e-3.
        rhow (int, float): water density, default rho is 1025 kg/m**3
           

    Synx & Desc:
        [UE,VE,wE] = ekman(lat, lon, u10, v10) estimates the zonal (UE, m**2/s) and meridional (VE, m**2/s)  
        Ekman layer transports along with vertical velocities (wE, m/s) associated with Ekman pumping. 
        Positive values of wE indicate upwelling. Inputs lat and lon must be 2D grids whose dimensions 
        zonal (u10, m/s) and meridional (v10, m/s) wind speeds taken 10 m above the surface. Wind speeds
        are automatically converted to wind stress before performing the Ekman transport calculations. 
        
        [UE,VE,wE] = ekman(lat, lon, u10, v10, Tau = True) estimates the zonal (UE, m**2/s) and meridional (VE, m**2/s)  
        Ekman layer transports along with vertical velocities (wE, m/s) associated with Ekman pumping. 
        Positive values of wE indicate upwelling. Inputs lat and lon must be 2D grids whose dimensions 
        zonal (u10, m/s) and meridional (v10, m/s) wind stress. 

        [UE,VE,wE] = ekman(...,Cd = Cd) specifies a drag coefficient for wind stress calculation. Default Cd is 1.25e-3. 

        [UE,VE,wE] = ekman(...,rhow = rhow) specifies water density. Default is 1025 kg/m**3.
    """
    # Calculate parameters for the input grid:
    # Grid dimensions:
    [dx, dy] = geoGrd.cdtdim(lat, lon)
    
    # Coriolis frequency:
    f = coriolisf(lat)
    
    if Tau:
        [TauX, TauY] = [u10, v10]
    
    else:
        # Estimate wind stress based on 10 m wind speed and drag coefficient:
        [TauX, TauY] = windstress(u10, v10, Cd)
    
    # Compute Ekman transports:
    # Zonal and meridional components:
    UE = TauY / (rhow * f)   # m**2/s
    VE = -TauX / (rhow * f)  # m**2/s
    
    # Compute vertical velocity:
    if UE.ndim == 3 and VE.ndim == 3:
        [_,_,UX] = np.gradient(UE)
        [_,VY,_] = np.gradient(VE)
    else:    
        [UX,_] = np.gradient(UE)
        [_,VY] = np.gradient(VE)
        
    UX = UX / dx
    VY = VY / dy
    wE = UX + VY
    
    return UE, VE, wE