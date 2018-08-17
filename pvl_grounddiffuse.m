function GR = pvl_grounddiffuse(SurfTilt,GHI,Albedo)
% PVL_GROUNDDIFFUSE Estimate diffuse irradiance from ground reflections given irradiance, albedo, and surface tilt 
%
% Syntax
%   GR = pvl_grounddiffuse(SurfTilt, GHI, Albedo)
%
% Description
%   Function to determine the portion of irradiance on a tilted surface due
%   to ground reflections. Any of the inputs may be vectors or scalars, 
%   as long as all vectors are of the same size.
%
% Inputs
%   SurfTilt - a scalar or vector of surface tilt angles in decimal degrees. 
%     If SurfTilt is a vector it must be of the same size as all other vector
%     inputs. SurfTilt must be >=0 and <=180. The tilt angle is defined as
%     degrees from horizontal (e.g. surface facing up = 0, surface facing
%     horizon = 90).
%   GHI - a scalar or vector of global horizontal irradiance in W/m^2. If GHI
%     is a vector it must be of the same size as all other vector inputs. 
%     GHI must be >=0.
%   Albedo - a scalar or vector for ground reflectance, typically 0.1-0.4 for
%     surfaces on Earth (land), may increase over snow, ice, etc. May also 
%     be known as the reflection coefficient. If Albedo is a vector it must
%     be of the same size as all other vector inputs. Must be >=0 and <=1.
%
% Outputs
%   GR is a column vector of ground reflected irradiances in W/m^2. 
%   The vector has the same number of elements as the input vector(s). 
%
% References
%   [1] Loutzenhiser P.G. et. al. "Empirical validation of models to compute
%   solar irradiance on inclined surfaces for building energy simulation"
%   2007, Solar Energy vol. 81. pp. 254-267
%
% See also PVL_DISC    PVL_PEREZ    PVL_REINDL1990    PVL_KLUCHER1979
% PVL_HAYDAVIES1980    PVL_ISOTROPICSKY    PVL_KINGDIFFUSE
%


p = inputParser;
p.addRequired('GHI',@(x) all(isnumeric(x) & isvector(x) & x>=0));
p.addRequired('Albedo', @(x) all(isnumeric(x) & isvector(x) & x>=0 & x<=1));
p.addRequired('SurfTilt', @(x) all(isvector(x) & isnumeric(x)));
p.parse(GHI,Albedo,SurfTilt);

GHI = GHI(:);
Albedo = Albedo(:);
SurfTilt = SurfTilt(:);

GR = GHI .* Albedo .* (1-cosd(SurfTilt)) .* 0.5;