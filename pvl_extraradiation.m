function Ea = pvl_extraradiation(doy)
% PVL_EXTRARADIATION Determine extraterrestrial radiation from day of year
%
% Syntax
%   Ea = pvl_extraradiation(doy)
%
% Description
%   Determine the amount of extraterrestrial solar radiation.
%
%   Output Ea is the extraterrestrial radiation present in watts per square meter
%   on a surface which is normal to the sun. Ea is of the same size as the
%   input doy.
%
%   Input doy is an array specifying the day of year. Valid values are >=1 and <367.
%
%
% Source
%   http://solardat.uoregon.edu/SolarRadiationBasics.html, Eqs. SR1 and SR2
%   SR1 	   	Partridge, G. W. and Platt, C. M. R. 1976. Radiative Processes in Meteorology and Climatology.
%   SR2 	   	Duffie, J. A. and Beckman, W. A. 1991. Solar Engineering of Thermal Processes, 2nd edn. J. Wiley and Sons, New York.
%
% See also PVL_DAYOFYEAR PVL_DISC

p =inputParser;
p.addRequired('doy',@(x) all(isnumeric(x) & x>=1 & x<367))
p.parse(doy)

B = 2*pi*doy/365;
Rfact2 = 1.00011 + 0.034221 .* cos(B)+ 0.00128.*sin(B)+ 0.000719.*cos(2*B)+0.000077.*sin(2*B);
Ea = 1367*Rfact2;
end