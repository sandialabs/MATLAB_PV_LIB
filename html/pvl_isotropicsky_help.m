%% pvl_isotropicsky 
% Determine sky diffuse irradiance on a tilted surface using the isotropic
% sky model.
%
%% Syntax
% |SkyDiffuse = pvl_isotropicsky(SurfTilt, DHI)|
%
%% Description
% The isotropic sky model [1], [2] regards the sky as a uniform source of diffuse
% irradiance. Thus the sky diffuse irradiance on a tilted surface can
% be found from the diffuse horizontal irradiance and the tilt angle of
% the surface.
%
%% Inputs 
%%
% * *|SurfTilt|* - a scalar or vector of surface tilt angles in decimal degrees. 
%     If SurfTilt is a vector it must be of the same size as all other vector
%     inputs. SurfTilt must be >=0 and <=180. The tilt angle is defined as
%     degrees from horizontal (e.g. surface facing up = 0, surface facing
%     horizon = 90)
% * *|DHI|* - a scalar or vector of diffuse horizontal irradiance in W/m^2. If DHI
%     is a vector it must be of the same size as all other vector inputs. 
%     DHI must be >=0.
%
%% Outputs
%%
% * *|SkyDiffuse|* - the diffuse component of the solar radiation  on an
%     arbitrarily tilted surface defined by the isotropic sky model as
%     given in Loutzenhiser et. al (2007) equation 3.
%     SkyDiffuse is the diffuse component ONLY and does not include the ground
%     reflected irradiance or the irradiance due to the beam.
%     SkyDiffuse is a column vector vector with a number of elements equal to
%     the input vector(s).
%
%% Example
% Calculate sky diffuse on a 30 deg tilted array when diffuse horizontal 
% irradiance equals 200 W/m^2.  
SkyDiffuse = pvl_isotropicsky(30, 200)
%% References
%%
% [1] Loutzenhiser P.G. et. al., 2007. Empirical validation of models to compute
%   solar irradiance on inclined surfaces for building energy simulation,
%   Solar Energy vol. 81. pp. 254-267.
%
% [2] Hottel, H.C., Woertz, B.B., 1942. Evaluation of flat-plate solar heat
%   collector. Trans. ASME 64, 91.
%
%% See Also 
% <pvl_perez_help.html |pvl_perez|>,
% <pvl_reindl1990_help.html |pvl_reindl1990|>, 
% <pvl_klucher1979_help.html |pvl_klucher1979|>, 
% <pvl_haydavies1980_help.html |pvl_haydavies1980|>, 
% <pvl_grounddiffuse_help.html |pvl_grounddiffuse|>, 
% <pvl_kingdiffuse_help.html |pvl_kingdiffuse|>
%   
%%
% Copyright 2014 Sandia National Laboratories

