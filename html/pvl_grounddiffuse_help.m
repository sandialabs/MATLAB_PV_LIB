%% pvl_grounddiffuse 
% Estimate diffuse irradiance from ground reflections given irradiance,
% albedo, and surface tilt.
%
%% Syntax
% |GR = pvl_grounddiffuse(SurfTilt, GHI, Albedo)|
%
%% Description
% Function to determine the portion of irradiance on a tilted surface due
% to ground reflections. Any of the inputs may be vectors or scalars, 
% as long as all vectors are of the same size.
%
%% Inputs
%%
% * *|SurfTilt|* - a scalar or vector of surface tilt angles in decimal degrees. 
%     If |SurfTilt| is a vector it must be of the same size as all other vector
%     inputs. |SurfTilt| must be >=0 and <=180. The tilt angle is defined as
%     degrees from horizontal (e.g. surface facing up = 0, surface facing
%     horizon = 90). 
% * *|GHI|* - a scalar or vector of global horizontal irradiance in W/m^2.
%     If |GHI| is a vector it must be of the same size as all other vector inputs. 
%     |GHI| must be >=0.
% * *|Albedo|* - a scalar or vector for ground reflectance, typically 0.1-0.4 for
%     surfaces on Earth (land), may increase over snow, ice, etc. May also 
%     be known as the reflection coefficient. If |Albedo| is a vector it must
%     be of the same size as all other vector inputs. Must be >=0 and <=1.
%
%% Outputs
%%
% * *|GR|* - a column vector of ground reflected irradiances in W/m^2. 
%   The vector has the same number of elements as the input vector(s). 
%
%% Example 1
% Calculate the ground reflected diffuse irradiance on an array tilted at 30 degrees
% from horizontal when GHI = 1,000 W/m^2 and albedo = 0.2.
GR = pvl_grounddiffuse(30, 1000, 0.2)
%% Example 2
% Calculate the ground reflected diffuse irradiance on an array tilted at 30 degrees
% from horizontal when GHI = 1,000 W/m^2 and albedo = 0.8 (typical of fresh
% snow).
GR = pvl_grounddiffuse(30, 1000, 0.8)
%% References
% [1] Loutzenhiser P.G. et al., 2007. Empirical validation of models to compute
%   solar irradiance on inclined surfaces for building energy simulation,
%   Solar Energy vol. 81. pp. 254-267.
%
%% See Also 
% <pvl_perez_help.html |pvl_perez|>,
% <pvl_reindl1990_help.html |pvl_reindl1990|>, 
% <pvl_klucher1979_help.html |pvl_klucher1979|>, 
% <pvl_haydavies1980_help.html |pvl_haydavies1980|>, 
% <pvl_isotropicsky_help.html |pvl_isotropicsky|>, 
% <pvl_kingdiffuse_help.html |pvl_kingdiffuse|>
%
%%
% Copyright 2014 Sandia National Laboratories


