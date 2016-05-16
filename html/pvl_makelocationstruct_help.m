%% pvl_makelocationstruct 
% Create a struct to define a site location.
%
%% Syntax
% |Location = pvl_makelocationstruct(latitude, longitude)|
%
% |Location = pvl_makelocationstruct(latitude, longitude, altitude)|
%
%% Description
% Creates a location struct for use with some PVLib functions. Site
% information includes latitude, longitude, and altitude (optional).
%
%% Inputs
% * *|latitude|* - latitude coordinates in decimal degrees. 
%             Latitude convention is positive north of the equator.
% * *|longitude|* - longitude coordinates in decimal degrees.
%           Longitude convention is positive east of the prime meridian.
% * *|altitude|* - (optional) surface elevation in meters above sea level.
%
% All inputs must be scalars.
%
%% Outputs
% * *|Location|* - Struct consisting of scalar components.
% * *|Location.latitude|* - Latitude in decimal degrees.
% * *|Location.longitude|* - Longitude in decimal degrees.
% * *|Location.altitude|* - Elevation in meters above sea level (if altitude provided). 
%
%% Example
%
latitude = 35.05;
longitude = -116.5;
altitude = 1200;
Location = pvl_makelocationstruct(latitude, longitude, altitude)
%% See Also 
% <pvl_ephemeris_help.html |pvl_ephemeris|>, 
% <pvl_maketimestruct_help.html |pvl_maketimestruct|>, 
% <pvl_alt2pres_help.html |pvl_alt2pres|>, 
% <pvl_pres2alt_help.html |pvl_pres2alt|>
%
%%
% Copyright 2014 Sandia National Laboratories