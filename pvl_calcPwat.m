function[Pwat] = pvl_calcPwat(T,RH)
% pvl_calcPwat Calculates precipitable water (cm) from ambient air temperature (C) and
%   relatively humidity (%)
%
% Syntax
%   [Pwat] = pvl_calcPwat(T,RH)
%
% Description
%   Calculates precipitable water (cm) from ambient air temperature (C) and
%   relatively humidity (%) using an empirical model [1]. The model was
%   developed by expanding Eq. 1 in [2]:
%           w = 0.1 H_v \rho_v 
%   using Eq. 2 in [2]
%           \rho_v = 216.7 RH/T e_s 
%   H_v is the apparant water vapor scale height (km). The expression for 
%   H_v is Eq. 4 in [2]:
%           H_v = 0.4976 + 1.5265.*T./273.15 + exp(13.6897.*T./273.15 - 14.9188.*(T./273.15).^3)
%   \rho_v is the surface water vapor density (g/m^3).  In the expression 
%   \rho_v,  e_s is the saturation water vapor pressure (millibar).  The
%   expression for e_s is Eq. 1 in [3]
%          e_s = exp(22.330 - 49.140.*(100./T) - 10.922.*(100./T).^2 - 0.39015.*T./100)
%
% Inputs:
%   T - ambient air temperature at the surface (C)
%   RH - relative humidity at the surface (%)
%
% Output: 
%   Pwat - precipitable water (cm)
%
% Reference:
% [1]   W. M. Keogh and A. W. Blakers, Accurate Measurement, Using Natural 
%           Sunlight, of Silicon Solar Cells, Prog. in Photovoltaics: Res. 
%           and Appl. 2004, vol 12, pp. 1-19 (DOI: 10.1002/pip.517)
% [2]   C. Gueymard, Analysis of Monthly Average Atmospheric Precipitable
%           Water and Turbidity in Canada and Northern United States,
%           Solar Energy vol 53(1), pp. 57-71, 1994.
% [3]   C. Gueymard, Assessment of the Accuracy and Computing Speed of
%           simplified saturation vapor equations using a new reference
%           dataset, J. of Applied Meteorology 1993, vol. 32(7), pp.
%           1294-1300.

T = T+273.15; %Convert to Kelvin

RH(RH>100 | RH<=0) = NaN; %Filter RH for unreasonable Values

% Eq. 1 from Keogh and Blakers
Pwat = 0.1.*(0.4976 + 1.5265.*T./273.15 + exp(13.6897.*T./273.15 - 14.9188.*(T./273.15).^3)).*...
    (216.7*RH./(100.*T).*exp(22.330 - 49.140.*(100./T) - 10.922.*(100./T).^2 - 0.39015.*T./100));

Pwat = max(Pwat,0.1); %Eliminates Pwat estimation less than 0.1cm
