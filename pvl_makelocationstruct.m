function Location = pvl_makelocationstruct(latitude, longitude, varargin)
% PVL_MAKELOCATIONSTRUCT Create a struct to define a site location
%
% Syntax
%   Location = pvl_makelocationstruct(latitude, longitude)
%   Location = pvl_makelocationstruct(latitude, longitude, altitude)
%
% Description
%   Creates a location struct for use with some PVLib functions. Site
%   information includes latitude, longitude, and optionally altitude.
%
%   Inputs latitude and longitude should be latitude and longitude
%   coordinates in decimal degrees. Latitude convention is positive north 
%   of the equator. Longitude convention is positive east of the prime meridian.
%   Input altitude is optional, and should be given in meters above sea
%   level.
%   All inputs must be scalars.
%
%   Output is a struct, Location, consisting of Location.latitude,
%   Location.longitude, and Location.altitude (if altitude provided). 
%   All output struct components are scalars.
%   
%
% See also PVL_EPHEMERIS PVL_MAKETIMESTRUCT PVL_ALT2PRES PVL_PRES2ALT
%
p = inputParser;
p.addRequired('latitude',@(x) all(isscalar(x) & isnumeric(x) & x<=90 & x>=-90));
p.addRequired('longitude', @(x) all(isscalar(x) & isnumeric(x) & x<=180 & x>=-180));
p.addOptional('altitude', 1.9321E100, @(x) all(isscalar(x) & isnumeric(x)));

p.parse(latitude,longitude,varargin{:});

Location.latitude = p.Results.latitude;
Location.longitude = p.Results.longitude;

defaultchecker = {'altitude'};
if ~any(strcmp(defaultchecker,p.UsingDefaults))
    Location.altitude = p.Results.altitude;
end
