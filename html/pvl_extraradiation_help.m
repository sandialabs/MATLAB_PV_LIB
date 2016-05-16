%% pvl_extraradiation 
% Determine extraterrestrial radiation from day of year
%
%% Syntax
% |Ea = pvl_extraradiation(doy)|
%
%% Description
% Determine the amount of extraterrestrial solar radiation.
%
%%  Inputs 
%%
% * *|doy|* - an array specifying the day of year. Valid values are >=1 and <367.
%
%%  Outputs 
%%
% * *|Ea|* - extraterrestrial radiation present in W/m^2
%  on a surface which is normal to the sun. |Ea| is of the same size as the
%  input |doy|.
%
%% Example
% Calculate extraterrestrial radiation for doy = 60;
Ea = pvl_extraradiation(60)
%% References
%%
% [1] http://solardat.uoregon.edu/SolarRadiationBasics.html, Eqs. SR1 from
% [2] and SR2 from [3].
%
% [2] Partridge, G. W. and Platt, C. M. R, Radiative Processes in
% Meteorology and Climatology, 1976.
%
% [3] Duffie, J. A. and Beckman, W. A., Solar Engineering of Thermal
% Processes, 2nd edn. J. Wiley and Sons, New York, 1991.
%
%% See Also 
% <pvl_date2doy_help |pvl_date2doy|>

%%
% Copyright 2014 Sandia National Laboratories
