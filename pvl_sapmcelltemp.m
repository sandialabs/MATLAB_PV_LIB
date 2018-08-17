function [Tcell Tmodule] = pvl_sapmcelltemp(E, E0, a, b, windspeed, Tamb, deltaT)
% PVL_SAPMCELLTEMP Estimate cell temperature from irradiance, windspeed, ambient temperature, and module parameters (SAPM)
%
% Syntax
%   Tcell = pvl_sapmcelltemp(E, E0, a, b, windspeed, Tamb, deltaT)
%   [Tcell Tmodule] = pvl_sapmcelltemp(E, E0, a, b, windspeed, Tamb, deltaT)
%
% Description
%   Estimate cell and module temperatures per the Sandia PV Array
%   Performance model (SAPM, SAND2004-3535), when given the incident
%   irradiance, wind speed, ambient temperature, and SAPM module
%   parameters.
%
% Inputs
%   E - Total incident irradiance in W/m^2. E must be a scalar or a vector
%     of the same size as windspeed, and Tamb. Must be >=0.
%   E0 - Reference irradiance used when determining delta T, in W/m^2. E0
%     must be a scalar. Must be >=0;
%   a - SAPM module parameter for establishing the upper limit for module 
%     temperature at low wind speeds and high solar irradiance (see SAPM
%     eqn. 11). Must be a scalar.
%   b - SAPM module parameter for establishing the rate at which the module
%     temperature drops as wind speed increases (see SAPM eqn. 11). Must be
%     a scalar.
%   windspeed - Wind speed in m/s at a height of 10 meters. windspeed must
%     be a scalar or a vector of the same size as E and Tamb. Must be >=0;
%   Tamb - Ambient dry bulb temperature in degrees C. Tamb must be a scalar
%     or a vector of the same size as windspeed and E. Must be >= -273.15.
%   deltaT - SAPM module parameter giving the temperature difference
%     between the cell and module back surface at the reference irradiance,
%     E0. Must be a numeric scalar >=0.
%  
% Outputs
%   Tcell - A column vector of cell temperatures in degrees C.
%   Tmodule - A column vector of module back temperature in degrees C.
%
% References
%   [1] King, D. et al, 2004, "Sandia Photovoltaic Array Performance Model", SAND Report
%   3535, Sandia National Laboratories, Albuquerque, NM
%
% See also PVL_SAPM
%
p = inputParser;
p.addRequired('E', @(x) all(isnumeric(x) & (x >= 0) & isvector(x)));
p.addRequired('E0', @(x) all(isnumeric(x) & isscalar(x) & (x >= 0)));
p.addRequired('a', @(x) all(isscalar(x) & isnumeric(x)));
p.addRequired('b', @(x) all(isscalar(x) & isnumeric(x)));
p.addRequired('windspeed', @(x) all(isvector(x) & isnumeric(x)));
p.addRequired('Tamb', @(x) all(isnumeric(x) & (x >= -273.15) & isvector(x)));
p.addRequired('deltaT', @(x) all(isnumeric(x) & (x >= 0) & isscalar(x)));
p.parse(E, E0, a, b, windspeed, Tamb, deltaT)

E = E(:);

if isscalar(p.Results.windspeed)
    windspeed =  p.Results.windspeed*ones(size(E));
else
    windspeed = p.Results.windspeed(:);
end

if isscalar(p.Results.Tamb)
    Tamb = p.Results.Tamb*ones(size(E));
else
    Tamb = p.Results.Tamb(:);
end

if (numel(E) ~= numel(windspeed)) || (numel(E) ~= numel(Tamb)) || (numel(E) ~= numel(windspeed))
    error(['Error in pvl_kingcelltemp. Inputs E, Tamb, and windspeed must '...
        'be scalars or vectors of same length.']);
end

Tmodule = E.*(exp(a+b.*windspeed))+Tamb;
Tcell = Tmodule + E./E0.*deltaT;