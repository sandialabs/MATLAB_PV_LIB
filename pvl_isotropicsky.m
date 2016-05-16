function SkyDiffuse = pvl_isotropicsky(SurfTilt, DHI)
% PVL_ISOTROPICSKY Determine diffuse irradiance from the sky on a tilted surface using isotropic sky model
%
% Syntax
%   SkyDiffuse = pvl_isotropicsky(SurfTilt, DHI)
%
% Description
%   Hottel and Woertz's model treats the sky as a uniform source of diffuse
%   irradiance. Thus the diffuse irradiance from the sky (ground reflected
%   irradiance is not included in this algorithm) on a tilted surface can
%   be found from the diffuse horizontal irradiance and the tilt angle of
%   the surface.
%
% Output:   
%   SkyDiffuse - the diffuse component of the solar radiation  on an
%     arbitrarily tilted surface defined by the isotropic sky model as
%     given in Loutzenhiser et. al (2007) equation 3.
%     SkyDiffuse is the diffuse component ONLY and does not include the ground
%     reflected irradiance or the irradiance due to the beam.
%     SkyDiffuse is a column vector vector with a number of elements equal to
%     the input vector(s).
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
%
% References
%   [1] Loutzenhiser P.G. et. al. "Empirical validation of models to compute
%   solar irradiance on inclined surfaces for building energy simulation"
%   2007, Solar Energy vol. 81. pp. 254-267
%   [2] Hottel, H.C., Woertz, B.B., 1942. Evaluation of flat-plate solar heat
%   collector. Trans. ASME 64, 91.
%
% See also    
%  PVL_REINDL1990  PVL_HAYDAVIES1980  PVL_PEREZ  PVL_KLUCHER1979
%  PVL_KINGDIFFUSE
%   
%

p = inputParser;
p.addRequired('SurfTilt', @(x) all(isnumeric(x) & x<=180 & x>=0 & isvector(x)));
p.addRequired('DHI', @(x) all(isnumeric(x) & isvector(x) & x>=0));
p.parse(SurfTilt, DHI);

SkyDiffuse = DHI *(1+ cosd(SurfTilt)) * 0.5;
SkyDiffuse = SkyDiffuse(:);