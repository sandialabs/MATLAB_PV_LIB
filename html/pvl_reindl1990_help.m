%% pvl_reindl1990 
% Determine sky diffuse irradiance on a tilted surface using
% Reindl's 1990 model.
%
%% Syntax
% |SkyDiffuse = pvl_reindl1990(SurfTilt, SurfAz, DHI, DNI, GHI, HExtra,
% SunZen, SunAz)|
%
%% Description
% Reindl's 1990 model ([1] equation 8, [2]) estimates the sky diffuse irradiance 
% on a tilted surface using the surface tilt angle, surface azimuth angle,
% diffuse horizontal irradiance, direct normal irradiance, global
% horizontal irradiance, extraterrestrial irradiance, sun zenith angle,
% and sun azimuth angle.
%

%% Inputs
%%
% * *|SurfTilt|* - a scalar or vector of surface tilt angles in decimal degrees. 
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
%     |DNI| is a vector it must be of the same size as all other vector inputs. 
%     |DNI| must be >=0.
% * *|GHI|* - a scalar or vector of global horizontal irradiance in W/m^2. 
%     If |GHI| is a vector it must be of the same size as all other vector inputs. 
%     |GHI| must be >=0.
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
%% Output
%%
% * *|SkyDiffuse|* - the sky diffuse irradiance on a tilted surface in W/m^2.
%     |SkyDiffuse| is a column vector vector with a number of elements equal to
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
SkyDiffuse = pvl_reindl1990(SurfTilt, SurfAz, DHI, DNI, GHI, HExtra,SunZen, SunAz)
%% References
%%
% [1] Loutzenhiser P.G. et. al. 2007. Empirical validation of models to compute
%   solar irradiance on inclined surfaces for building energy simulation,
%   Solar Energy vol. 81. pp. 254-267.
%
% [2] Reindl, D.T., Beckmann, W.A., Duffie, J.A., 1990. Evaluation of hourly
%   tilted surface radiation models. Solar Energy 45 (1), 9–17.
%
%% See Also 
% <pvl_perez_help.html |pvl_perez|>,
% <pvl_grounddiffuse_help.html |pvl_grounddiffuse|>, 
% <pvl_klucher1979_help.html |pvl_klucher1979|>, 
% <pvl_haydavies1980_help.html |pvl_haydavies1980|>, 
% <pvl_isotropicsky_help.html |pvl_isotropicsky|>, 
% <pvl_kingdiffuse_help.html |pvl_kingdiffuse|>,
% <pvl_extraradiation_help.html |pvl_extraradiation|>
%
%%
% Copyright 2014 Sandia National Laboratories