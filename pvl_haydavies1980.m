function SkyDiffuse = pvl_haydavies1980(SurfTilt, SurfAz, DHI, DNI, HExtra, SunZen, SunAz)
% PVL_HAYDAVIES1980 Determine diffuse irradiance from the sky on a tilted surface using Hay & Davies' 1980 model
%
% Syntax
%   SkyDiffuse = pvl_haydavies1980(SurfTilt, SurfAz, DHI, DNI, HExtra, SunZen, SunAz)
%
% Description
%   Hay and Davies' 1980 model determines the diffuse irradiance from the sky
%   (ground reflected irradiance is not included in this algorithm) on a
%   tilted surface using the surface tilt angle, surface azimuth angle,
%   diffuse horizontal irradiance, direct normal irradiance, 
%   extraterrestrial irradiance, sun zenith angle, and sun azimuth angle.
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
%     DHI must be >=0.
%   DNI - a scalar or vector of direct normal irradiance in W/m^2. If DNI
%     is a vector it must be of the same size as all other vector inputs. 
%     DNI must be >=0.
%   HExtra - a scalar or vector of extraterrestrial normal irradiance in 
%     W/m^2. If HExtra is a vector it must be of the same size as 
%     all other vector inputs. HExtra must be >=0.
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
%     arbitrarily tilted surface defined by the Hay & Davies model as given in
%     Loutzenhiser et. al (2007) equation 7.
%     SkyDiffuse is the diffuse component ONLY and does not include the ground
%     reflected irradiance or the irradiance due to the beam.
%     SkyDiffuse is a column vector vector with a number of elements equal to
%     the input vector(s).
%
% References
%   [1] Loutzenhiser P.G. et. al. "Empirical validation of models to compute
%   solar irradiance on inclined surfaces for building energy simulation"
%   2007, Solar Energy vol. 81. pp. 254-267
%   [2] Hay, J.E., Davies, J.A., 1980. Calculations of the solar radiation incident
%   on an inclined surface. In: Hay, J.E., Won, T.K. (Eds.), Proc. of First
%   Canadian Solar Radiation Data Workshop, 59. Ministry of Supply
%   and Services, Canada.
%
% See also PVL_EPHEMERIS   PVL_EXTRARADIATION   PVL_ISOTROPICSKY
%       PVL_REINDL1990   PVL_PEREZ PVL_KLUCHER1979   PVL_KINGDIFFUSE
%        PVL_SPA
%
%
%
%

p=inputParser;
p.addRequired('SurfTilt', @(x) (isnumeric(x) && all(x<=180) && all(x>=0) && isvector(x)));
p.addRequired('SurfAz', @(x) isnumeric(x) && all(x<=360) && all(x>=0) && isvector(x));
p.addRequired('DHI', @(x) (isnumeric(x) && isvector(x) && all((x>=0) | isnan(x))));
p.addRequired('DNI', @(x) isnumeric(x) && isvector(x) && all((x>=0) | isnan(x)));
p.addRequired('HExtra', @(x) isnumeric(x) && isvector(x) && all((x>=0) | isnan(x)));
p.addRequired('SunZen', @(x) isnumeric(x) && all(x<=180) && all((x>=0) | isnan(x)) && isvector(x));
p.addRequired('SunAz', @(x) (isnumeric(x) && all(x<=360) && all((x>=0) | isnan(x)) && isvector(x)));
p.parse(SurfTilt, SurfAz, DHI, DNI, HExtra, SunZen, SunAz);

%COSTT is the cosine of the angle of incidence between beam of the sun
%(defined by SunAz and SunZen)and the normal vector to a surface (defined
% by SurfAz and SurfTilt). 
%
% Dec 2012: A bug was identified by Rob Andrews (Queens University) in this equation in PV_LIB
% Version 1.0.  Fixed in Version 1.1.
COSTT = cosd(SurfTilt).*cosd(SunZen) + sind(SurfTilt).* ...
    sind(SunZen).*cosd(SunAz-SurfAz); 
% RB is the ratio of the beam irradiance on the tilted surface to the beam
% irradiance on a horizontal surface, the denominator is set to a minimum
% of 0.01745 corresponding to a zenith angle of 89 degrees
RB = max(COSTT,0)./max(cosd(SunZen),0.01745);
% AI is the Anisotropy index which represents the transmittance through
% atmosphere for beam radiation
AI = DNI./HExtra;

% This is the POAskydiffuse calculation I generated from the Loutzenhiser et al.
% (2007) paper, equation 7. Note that I have removed the beam and ground
% reflectance portion of the equation and this generates ONLY the diffuse
% radiation from the sky and circumsolar, so the form of the equation
% varies slightly from equation 7.
SkyDiffuse = DHI .* (AI .* RB + (1-AI) .* 0.5 .* (1 + cosd(SurfTilt)));
SkyDiffuse = SkyDiffuse(:);

end