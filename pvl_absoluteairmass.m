function AMa = pvl_absoluteairmass(AMrelative,pressure)
% PVL_ABSOLUTEAIRMASS Determine absolute (pressure corrected) airmass from relative airmass and pressure
%
% Syntax
% AMa = pvl_absoluteairmass(AMrelative, pressure)
%
% Description
% 
%   Gives the airmass for locations not at sea-level (i.e. not at standard
%   pressure). The input argument "AMrelative" is the relative airmass. The
%   input argument "pressure" is the pressure (in Pascals) at the location
%   of interest and must be greater than 0. The calculation for
%   absolute airmass is:
%   absolute airmass = (relative airmass)*pressure/101325
%
% Inputs:   
%   AMrelative - The airmass at sea-level.  This can be calculated using the 
%     PV_LIB function pvl_relativeairmass. 
%   pressure - a scalar or vector of values providing the site pressure in
%     Pascal. If pressure is a vector it must be of the same size as all
%     other vector inputs. pressure must be >=0. Pressure may be measured
%     or an average pressure may be calculated from site altitude.
%
% Output:   
%   AMa - Absolute (pressure corrected) airmass
%   
% References
%   [1] C. Gueymard, "Critical analysis and performance assessment of 
%   clear sky solar irradiance models using theoretical and measured data,"
%   Solar Energy, vol. 51, pp. 121-138, 1993.
%
% See also PVL_RELATIVEAIRMASS

p=inputParser;
p.addRequired('AMrelative', @(x) all(isnumeric(x) | isnan(x)));
p.addRequired('pressure', @(x) all(isnumeric(x) & x>=0));
p.parse(AMrelative,pressure);

AMa = AMrelative.*pressure/101325;