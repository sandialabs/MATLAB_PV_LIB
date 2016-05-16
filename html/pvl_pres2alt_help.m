%% pvl_pres2alt 
% Determine altitude (m) from site pressure (Pa).
%
%% Syntax
% |altitude = pvl_pres2alt(pressure)| 
%
%% Description
% |altitide = pvl_pres2alt(pressure)| determines the altitude (in meters above sea level) of a 
%   site on Earth's surface given its atmospheric pressure (in Pascals). 
%
%% Assumptions
%
% * Pressure at mean sea level = 101325 Pa.
% * Temperature at zero altitude = 288.15 K.
% * Gravitational acceleration = 9.80665 m/s^2.
% * Lapse rate = -6.5E-3 K/m.
% * Gas constant for air = 287.053 J/(kg*K).
% * Relative Humidity = 0%.
%
%% Inputs
%%
% * *|pressure|* - atmospheric pressure (in Pascals) of a 
% site on Earth's surface given its altitude
%
%% Outputs
%%
% * *|altitude|* - altitude (in meters above sea level) as a vector of the same size 
%   as |pressure|.
%
%% Example
Alt1 = pvl_pres2alt(101325) % Sea level, result should be near zero 
Alt2 = pvl_pres2alt(80000) % Approximate pressure near 2,000 meters abobe sea level
%% References
% Portland State Aerospace Society, 2004. A Quick Derivation relating 
% altitude to air pressure, Version 1.03,
% <http://psas.pdx.edu/RocketScience/PressureAltitude_Derived.pdf>.
%
%% See Also
% <pvl_alt2pres_help.html |pvl_alt2pres|>,    
% <pvl_makelocationstruct_help.html |pvl_makelocationstruct|>
%%
% Copyright 2014 Sandia National Laboratories
