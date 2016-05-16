function ClearSkyGHI = pvl_clearsky_haurwitz(ApparentZenith)
% PVL_CLEARSKY_HAURWITZ Determine clear sky GHI from Haurwitz model
%
% Syntax
%   [ClearSkyGHI]= pvl_clearsky_haurwitz(ApparentZenith)
%   
% Description
%   Implements the Haurwitz clear sky model for global horizontal
%   irradiance (GHI) as presented in [1, 2]. A report on clear
%   sky models found the Haurwitz model to have the best performance of
%   models which require only zenith angle [3].
%
% Input Parameters:
%   ApparentZenith - The apparent (refraction corrected) sun zenith angle
%       in degrees.
%
% Output:   
%   ClearSkyGHI - the modeled global horizonal irradiance in W/m^2 provided
%      by the Haurwitz clear-sky model.
%
%   Initial implementation of this algorithm by Matthew Reno.
%
% Sources:
%
% [1] B. Haurwitz, "Insolation in Relation to Cloudiness and Cloud 
%     Density," Journal of Meteorology, vol. 2, pp. 154-166, 1945.
%
% [2] B. Haurwitz, "Insolation in Relation to Cloud Type," Journal of 
%     Meteorology, vol. 3, pp. 123-124, 1946.
%
% [3] M. Reno, C. Hansen, and J. Stein, "Global Horizontal Irradiance Clear
%     Sky Models: Implementation and Analysis", Sandia National
%     Laboratories, SAND2012-2389, 2012.
%
% See also
%   PVL_MAKETIMESTRUCT    PVL_MAKELOCATIONSTRUCT   PVL_EPHEMERIS   PVL_SPA
%   PVL_CLEARSKY_INEICHEN

p = inputParser;
p.addRequired('ApparentZenith', @(x) (all(isnumeric(x) & x<=180 & x>=0 & isvector(x))));
p.parse(ApparentZenith);

% Haurwitz uses 94.4 cal/cm^2/hr, this converts to approximately 1098 W/m^2
ClearSkyGHI = 1098.* cosd(ApparentZenith) .* exp(-0.059 ./ cosd(ApparentZenith));
ClearSkyGHI(ClearSkyGHI<0)=0;



end