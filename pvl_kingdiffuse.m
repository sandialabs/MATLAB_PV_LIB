function SkyDiffuse = pvl_kingdiffuse(SurfTilt, DHI, GHI, SunZen)
% PVL_KINGDIFFUSE Determine the total diffuse irradiance on a tilted surface using
% the King model.
%
% Syntax
%   POADiffuse = pvl_kingdiffuse(SurfTilt, DHI, GHI, SunZen)
%
% Description
%   The Kingis model (developed by David L. King at Sandia National Laboratories)
%   determines the total diffuse irradiance (sky diffuse and ground reflected irradiance) 
%   on a tilted surface using the surface tilt angle, diffuse horizontal
%   irradiance, global horizontal irradiance, and sun zenith angle. 
%
% Inputs:   
%   SurfTilt - a scalar or vector of surface tilt angles in decimal degrees. 
%     If SurfTilt is a vector it must be of the same size as all other vector
%     inputs. SurfTilt must be >=0 and <=180. The tilt angle is defined as
%     degrees from horizontal (e.g. surface facing up = 0, surface facing
%     horizon = 90)
%   DHI - a scalar or vector of diffuse horizontal irradiance in W/m^2. If DHI
%     is a vector it must be of the same size as all other vector inputs. 
%     DHI must be >=0.
%   GHI - a scalar or vector of global horizontal irradiance in W/m^2. If GHI
%     is a vector it must be of the same size as all other vector inputs. 
%     GHI must be >=0.
%   SunZen - a scalar or vector of apparent (refraction-corrected) zenith
%     angles in decimal degrees. If SunZen is a vector it must be of the
%     same size as all other vector inputs. SunZen must be >=0 and <=180.
%
% Output:   
%   POADiffuse - the total diffuse irradiance on a tilted surface.
%
% References
%   [1] Lave, M., Hayes, W., Pohl, A., Hansen, C., 2015. Evaluation of
%   Global Horizontal Irradiance to Plane-of-Array Irradiance Models at
%   Locations Across the United States, J. of Photovoltaics v5(2), pp.
%   597-606.
%
% See also PVL_EPHEMERIS   PVL_EXTRARADIATION   PVL_ISOTROPICSKY
%       PVL_HAYDAVIES1980   PVL_PEREZ PVL_KLUCHER1979   PVL_REINDL1990
%
p = inputParser;
p.addRequired('SurfTilt', @(x) all(isnumeric(x) & x<=180 & x>=0 & isvector(x)));
p.addRequired('DHI', @(x) all(isnumeric(x) & isvector(x) & x>=0));
p.addRequired('GHI', @(x) all(isnumeric(x) & isvector(x) & x>=0));
p.addRequired('SunZen', @(x) all(isnumeric(x) & x<=180 & x>=0 & isvector(x)));
p.parse(SurfTilt, DHI, GHI, SunZen);

SkyDiffuse = DHI.*(1+cosd(SurfTilt))./2 +...
    GHI.*(.012*SunZen-.04).*(1-cosd(SurfTilt))./2;
SkyDiffuse = SkyDiffuse(:);