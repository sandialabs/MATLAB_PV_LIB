function SkyDiffuse = pvl_klucher1979(SurfTilt, SurfAz, DHI, GHI, SunZen, SunAz)
% PVL_KLUCHER1979 Determine diffuse irradiance from the sky on a tilted surface using Klucher's 1979 model
%
% Syntax
%   SkyDiffuse =  pvl_klucher1979(SurfTilt, SurfAz, DHI, GHI, SunZen, SunAz)
%
% Description
%   Klucher's 1979 model determines the diffuse irradiance from the sky
%   (ground reflected irradiance is not included in this algorithm) on a
%   tilted surface using the surface tilt angle, surface azimuth angle,
%   diffuse horizontal irradiance, direct normal irradiance, global
%   horizontal irradiance, extraterrestrial irradiance, sun zenith angle,
%   and sun azimuth angle.
%
% Inputs:   
%   SurfTilt - a scalar or vector of surface tilt angles in decimal degrees. 
%     If SurfTilt is a vector it must be of the same size as all other vector
%     inputs. SurfTilt must be >=0 and <=180. The tilt angle is defined as
%     degrees from horizontal (e.g. surface facing up = 0, surface facing
%     horizon = 90)
%   SurfAz - a scalar or vector of surface azimuth angles in decimal degrees. 
%     If SurfAz is a vector it must be of the same size as all other vector
%     inputs. SurfAz must be >=0 and <=360. The Azimuth convention is defined
%     as degrees east of north (e.g. North = 0, East = 90, West = 270).
%   DHI - a scalar or vector of diffuse horizontal irradiance in W/m^2. If DHI
%     is a vector it must be of the same size as all other vector inputs. 
%     DHI must be >=0. Values of DHI which are < GHI will be set to GHI.
%   GHI - a scalar or vector of global horizontal irradiance in W/m^2. If GHI
%     is a vector it must be of the same size as all other vector inputs. 
%     GHI must be >=0.
%   SunZen - a scalar or vector of apparent (refraction-corrected) zenith
%     angles in decimal degrees. If SunZen is a vector it must be of the
%     same size as all other vector inputs. SunZen must be >=0 and <=180.
%   SunAz - a scalar or vector of sun azimuth angles in decimal degrees. 
%     If SunAz is a vector it must be of the same size as all other vector
%     inputs. SunAz must be >=0 and <=360. The Azimuth convention is defined
%     as degrees east of north (e.g. North = 0, East = 90, West = 270).
%
% Output:   
%   SkyDiffuse - the diffuse component of the solar radiation  on an
%     arbitrarily tilted surface defined by the Klucher model as given in
%     Loutzenhiser et. al (2007) equation 4.
%     SkyDiffuse is the diffuse component ONLY and does not include the ground
%     reflected irradiance or the irradiance due to the beam.
%     SkyDiffuse is a column vector vector with a number of elements equal to
%     the input vector(s).
%
% References
%   [1] Loutzenhiser P.G. et. al. "Empirical validation of models to compute
%   solar irradiance on inclined surfaces for building energy simulation"
%   2007, Solar Energy vol. 81. pp. 254-267
%   [2] Klucher, T.M., 1979. Evaluation of models to predict insolation on tilted
%   surfaces. Solar Energy 23 (2), 111–114.
%
% See also PVL_EPHEMERIS   PVL_EXTRARADIATION   PVL_ISOTROPICSKY
%       PVL_HAYDAVIES1980   PVL_PEREZ  PVL_REINDL1990   PVL_KINGDIFFUSE
%
%
%
%

p = inputParser;
p.addRequired('SurfTilt', @(x) all(isnumeric(x) & x<=180 & x>=0 & isvector(x)));
p.addRequired('SurfAz', @(x) all(isnumeric(x) & x<=360 & x>=0 & isvector(x)));
p.addRequired('DHI', @(x) all(isnumeric(x) & isvector(x) & x>=0));
p.addRequired('GHI', @(x) all(isnumeric(x) & isvector(x) & x>=0));
p.addRequired('SunZen', @(x) all(isnumeric(x) & x<=180 & x>=0 & isvector(x)));
p.addRequired('SunAz', @(x) all(isnumeric(x) & x<=360 & x>=0 & isvector(x)));
p.parse(SurfTilt, SurfAz, DHI, GHI, SunZen, SunAz);

GHI(GHI<DHI) = DHI(GHI<DHI);    %Global horizontal should be >= Diffuse horizontal
GHI(GHI < 1E-6) = 1E-6; % Prevent division by 0 
% Dec 2012: A bug was identified by Rob Andrews (Queens University) in this equation in PV_LIB
% Version 1.0.  Fixed in Version 1.1.
COSTT = cosd(SurfTilt).*cosd(SunZen) + sind(SurfTilt).* ...
    sind(SunZen).*cosd(SunAz-SurfAz);

F = 1 - ((DHI./GHI).^2);
SkyDiffuse = DHI ...
    .* (0.5 .* (1 + cosd(SurfTilt))) ...
    .* (1 + F .* ((sind(SurfTilt ./ 2)).^3)) ...
    .* (1 + F .* ((COSTT).^2) .* ((sind(SunZen)).^3));

SkyDiffuse = SkyDiffuse(:); % Make the column vector, regardless of input vector type

