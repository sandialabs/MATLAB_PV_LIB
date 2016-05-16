function SkyDiffuse = pvl_reindl1990(SurfTilt, SurfAz, DHI, DNI, GHI, HExtra, SunZen, SunAz)
% PVL_REINDL1990 Determine diffuse irradiance from the sky on a tilted surface using Reindl's 1990 model
%
% Syntax
%   SkyDiffuse = pvl_reindl1990(SurfTilt, SurfAz, DHI, DNI, GHI, HExtra, SunZen, SunAz)
%
% Description
%   Reindl's 1990 model determines the diffuse irradiance from the sky
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
%     DHI must be >=0.
%   DNI - a scalar or vector of direct normal irradiance in W/m^2. If DNI
%     is a vector it must be of the same size as all other vector inputs. 
%     DNI must be >=0.
%   GHI - a scalar or vector of global horizontal irradiance in W/m^2. If GHI
%     is a vector it must be of the same size as all other vector inputs. 
%     GHI must be >=0.
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
%     arbitrarily tilted surface defined by the Reindl model as given in
%     Loutzenhiser et. al (2007) equation 8.
%     SkyDiffuse is the diffuse component ONLY and does not include the ground
%     reflected irradiance or the irradiance due to the beam.
%     SkyDiffuse is a column vector vector with a number of elements equal to
%     the input vector(s).
%
% References
%   [1] Loutzenhiser P.G. et. al. "Empirical validation of models to compute
%   solar irradiance on inclined surfaces for building energy simulation"
%   2007, Solar Energy vol. 81. pp. 254-267
%   [2] Reindl, D.T., Beckmann, W.A., Duffie, J.A., 1990a. Diffuse fraction
%   correlations. Solar Energy 45 (1), 1–7.
%   [3] Reindl, D.T., Beckmann, W.A., Duffie, J.A., 1990b. Evaluation of hourly
%   tilted surface radiation models. Solar Energy 45 (1), 9–17.
%
% See also PVL_EPHEMERIS   PVL_EXTRARADIATION   PVL_ISOTROPICSKY
%       PVL_HAYDAVIES1980   PVL_PEREZ PVL_KLUCHER1979   PVL_KINGDIFFUSE
%
%
%
%

p = inputParser;
p.addRequired('SurfTilt', @(x) all(isnumeric(x) & x<=180 & x>=0 & isvector(x)));
p.addRequired('SurfAz', @(x) all(isnumeric(x) & x<=360 & x>=0 & isvector(x)));
p.addRequired('DHI', @(x) all(isnumeric(x) & isvector(x) & x>=0));
p.addRequired('DNI', @(x) all(isnumeric(x) & isvector(x) & x>=0));
p.addRequired('GHI', @(x) all(isnumeric(x) & isvector(x) & x>=0));
p.addRequired('HExtra', @(x) all(isnumeric(x) & isvector(x) & x>=0));
p.addRequired('SunZen', @(x) all(isnumeric(x) & x<=180 & x>=0 & isvector(x)));
p.addRequired('SunAz', @(x) all(isnumeric(x) & x<=360 & x>=0 & isvector(x)));
p.parse(SurfTilt, SurfAz, DHI, DNI, GHI, HExtra, SunZen, SunAz);



% Function to calculate Reindl Model for Diffuse POA
small = 0.000001; %from TRANSYS, just used to make sure that we don't divide by 0
% Dec 2012: A bug was identified by Rob Andrews (Queens University) in this equation in PV_LIB
% Version 1.0.  Fixed in Version 1.1.
COSTT = cosd(SurfTilt).*cosd(SunZen) + sind(SurfTilt).*sind(SunZen).*...
    cosd(SunAz-SurfAz);
RB = max(COSTT,0)./max(cosd(SunZen),0.01745);
AI = DNI./HExtra;

GHI(GHI<small)=small; % prevent division by zero

%See Reindl's "Evaluation of Hourly Tilted Surface Radiation Models" (1990)
%for more information on the modulation factor, F.
HB = DNI.*cosd(SunZen);
HB(HB<0) = 0; % Don't take the sqaure root of a negative number
GHI(GHI<0) = 0; % Don't take the sqaure root of a negative number
F = sqrt(HB./GHI);
SCUBE = (sind(SurfTilt.*0.5)).^3;

% This is the POAskydiffuse calculation I generated from the Loutzenhiser et al.
% (2007) paper, equation 8. Note that I have removed the beam and ground
% reflectance portion of the equation and this generates ONLY the diffuse
% radiation from the sky and circumsolar, so the form of the equation
% varies slightly from equation 8.
SkyDiffuse = DHI .* (AI .* RB + (1-AI) .* 0.5 .* (1 + cosd(SurfTilt)) .* ...
    (1 + F .* SCUBE));
SkyDiffuse = SkyDiffuse(:); % Make the column vector, regardless of input vector type

