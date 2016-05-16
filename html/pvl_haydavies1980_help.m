%% pvl_haydavies1980 
% Determine diffuse irradiance from the sky on a tilted surface using Hay &
% Davies' 1980 model.
%
%% Syntax
% |SkyDiffuse = pvl_haydavies1980(SurfTilt, SurfAz, DHI, DNI, HExtra,
% SunZen, SunAz)|
%
%% Description
% Hay and Davies' 1980 model determines the sky diffuse irradiance on a
% tilted surface using the surface tilt angle, surface azimuth angle,
% diffuse horizontal irradiance, direct normal irradiance, 
% extraterrestrial irradiance, sun zenith angle, and sun azimuth angle. 
%
%% Inputs
%%
% * *|SurfTilt|*- a scalar or vector of surface tilt angles in decimal degrees. 
%     If |SurfTilt| is a vector it must be of the same size as all other vector
%     inputs. |SurfTilt| must be >=0 and <=180. The tilt angle is defined as
%     degrees from horizontal (e.g. surface facing up = 0, surface facing
%     horizon = 90)
% * *|SurfAz|* - a scalar or vector of surface azimuth angles in decimal degrees. 
%     If |SurfAz| is a vector it must be of the same size as all other vector
%     inputs. |SurfAz| must be >=0 and <=360. The azimuth convention is defined
%     as degrees east of north (e.g. North = 0, East = 90, West = 270).
% * *|DHI|* - a scalar or vector of diffuse horizontal irradiance in W/m^2.
%     If |DHI| is a vector it must be of the same size as all other vector inputs. 
%     |DHI| must be >=0.
% * *|DNI|* - a scalar or vector of direct normal irradiance in W/m^2. If
%    |DNI| is a vector it must be of the same size as all other vector inputs. 
%    |DNI| must be >=0.
% * *|HExtra|* - a scalar or vector of extraterrestrial normal irradiance in 
%     W/m^2. If |HExtra| is a vector it must be of the same size as 
%     all other vector inputs. |HExtra| must be >=0.
% * *|SunZen|* - a scalar or vector of apparent (refraction-corrected) zenith
%     angles in decimal degrees. If |SunZen| is a vector it must be of the
%     same size as all other vector inputs. |SunZen| must be >=0 and <=180.
% * *|SunAz|* - a scalar or vector of sun azimuth angles in decimal degrees. 
%     If |SunAz| is a vector it must be of the same size as all other vector
%     inputs. |SunAz| must be >=0 and <=360. The azimuth convention is defined
%     as degrees east of north (e.g. North = 0, East = 90, West = 270).
%
%% Outputs
% * *|SkyDiffuse|* - the diffuse component of the solar radiation  on an
%     arbitrarily tilted surface defined by the Hay & Davies model [2] as 
%     summarized in [1] Equation 7.
%     SkyDiffuse is the sky diffuse component ONLY and does not include the ground
%     reflected diffuse irradiance or direct irradiance.
%     SkyDiffuse is a column vector vector with a number of elements equal to
%     the input vector(s).
%
%% Example
%
SurfTilt = 30;
SurfAz = 180;
DHI = 47;
DNI = 969;
GHI = 473;
HExtra = pvl_extraradiation(60);
SunZen = 60;
SunAz = 161;
pvl_haydavies1980(SurfTilt, SurfAz, DHI, DNI, HExtra, SunZen, SunAz)

%% References
%%
% [1] Loutzenhiser P.G. et. al., 2007. Empirical validation of models to compute
%   solar irradiance on inclined surfaces for building energy simulation,
%   Solar Energy vol. 81. pp. 254-267
%
% [2] Hay, J.E., Davies, J.A., 1980. Calculations of the solar radiation incident
%   on an inclined surface. In: Hay, J.E., Won, T.K. (Eds.), Proc. of First
%   Canadian Solar Radiation Data Workshop, 59. Ministry of Supply
%   and Services, Canada.
%
%% See Also 
% <pvl_perez_help.html |pvl_perez|>,
% <pvl_reindl1990_help.html |pvl_reindl1990|>, 
% <pvl_kingdiffuse_help.html |pvl_kingdiffuse|>, 
% <pvl_klucher1979_help.html |pvl_klucher1979|>, 
% <pvl_grounddiffuse_help.html |pvl_grounddiffuse|>, 
% <pvl_ephemeris_help.html |pvl_ephemeris|>   
%
%%
% Copyright 2014 Sandia National Laboratories

