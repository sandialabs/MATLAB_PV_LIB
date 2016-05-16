%% pvl_absoluteairmass
% Determine absolute (pressure corrected) airmass from relative airmass and pressure
%
%% Syntax
% |AMa = pvl_absoluteairmass(AMrelative, pressure)|
%
%% Description
% Gives the airmass for locations not at sea-level (i.e. not at standard
% pressure). The input argument |AMrelative| is the relative airmass. The
% input argument |pressure| is the pressure (in Pascals) at the location
% of interest and must be greater than 0. The calculation for
% absolute airmass is:
%
% $$AMa = AM \times {\frac{Pressure}{101325}}$$
%
%% Inputs 
%%
% * *|AMrelative|* - The airmass at sea-level.  This can be calculated using the 
%     PV_LIB function <pvl_relativeairmass_help.html |pvl_relativeairmass|>. 
% * *|pressure|* - a scalar or vector of values providing the site pressure in
%     Pascal. If pressure is a vector it must be of the same size as all
%     other vector inputs. pressure must be >=0. Pressure may be measured
%     or an average pressure may be calculated from site altitude.
%
%% Outputs 
%%
% * *|AMa|* - Absolute (pressure corrected) airmass
%% Example
% Adjust airmass for elevation 
AM = 1;
pressure = 80000; %Approximate air pressure at 2000 meters above sea level
AMa = pvl_absoluteairmass(AM, pressure)

%% References
% [1] C. Gueymard, "Critical analysis and performance assessment of 
% clear sky solar irradiance models using theoretical and measured data,"
% Solar Energy, vol. 51, pp. 121-138, 1993.
%
%% See Also 
% <pvl_relativeairmass_help.html |pvl_relativeairmass|>
%
%%
% Copyright 2014 Sandia National Laboratories

