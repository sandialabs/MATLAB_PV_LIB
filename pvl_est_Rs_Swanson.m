function [Rs, errest] = pvl_est_Rs_Swanson(V1, I1, V2, I2, delI, varargin)
% PVL_EST_RS_SWANSON estimates Rs using the method of Swanson (as described
% in Pysch [1]).
%
% Syntax
% [Rs, errest] = pvl_est_Rs_Swanson(V1, I1, V2, I2, delI)
% [Rs, errest] = pvl_est_Rs_Swanson(V1, I1, V2, I2, delI, Rsh1, Rsh2, nNsVth, Io)
%
% Description
%   The method of Swanson (see [1]) estimates series resistance using 
%   points from two IV curves at unequal but similar irradiance, e.g.,
%   at 950 W/m2 and 1000 W/m2. The cell temperature is assumed to be the 
%   same for both IV curves. A point is selected on each IV curve where 
%   the current is less than Isc by the input delI. If optional arguments are
%   provided, the difference between the returned Rs value and the Rs
%   parameter for the single diode equation is estimated, see [2].
%
% Inputs:
%   V1, I1 - vectors of voltage and current for the first IV curve.
%   V2, I2 - vectors of voltage and current for the first IV curve.
%   delI - offset from Isc to use
%   Rsh1 - (optional) shunt resistance, in ohms, for the first IV curve.
%   Rsh2 - (optional) shunt resistance, in ohms, for the second IV curve.
%   nNsVth - (optional) the product n (diode factor) x Ns (cells in series)
%     x Vth (thermal voltage per cell) for both IV curves.
%   Io - (optional) the dark current, in A, for both IV curves.
%
% Output:
%   Rs - the series resistance value in ohms.
%   errest - the estimated difference between Rs and the series resistance
%     parameter for the single diode equation.
%
% References
%   [1] D. Pysch, A. Mette, S. W. Glunz, “A review and comparison of 
%   different methods to determine the series resistance of solar cells",
%   Solar. Energy Materials and Cells 91, pp. 1698-1706, 2007.
%
%   [2] C. Hansen and B. King, "Determining series resistance for
%   equivalent circuit models of a PV module", in 45th IEEE Photovoltaic
%   Specialist Conference, Waikoloa, HI, 2018.

p = inputParser;
addRequired(p, 'V1', @isvector);
addRequired(p, 'I1', @isvector);
addRequired(p, 'V2', @isvector);
addRequired(p, 'I2', @isvector);
addRequired(p, 'delI', @isnumeric);

addOptional(p, 'Rsh1', NaN, @isnumeric);
addOptional(p, 'Rsh2', NaN, @isnumeric);
addOptional(p, 'nNsVth', NaN, @isnumeric);
addOptional(p, 'Io', NaN, @isnumeric);

parse(p, V1, I1, V2, I2, delI, varargin{:});
V1 = p.Results.V1;
I1 = p.Results.I1;
V2 = p.Results.V2;
I2 = p.Results.I2;
delI = p.Results.delI;

Rsh1 = p.Results.Rsh1;
Rsh2 = p.Results.Rsh2;
nNsVth = p.Results.nNsVth;
Io = p.Results.Io;

% calculate current on each IV curve
pI1 = I1(1) - delI;
pI2 = I2(1) - delI;

if pI1>0 && pI2>0
    
    % find voltage on each IV curve
    pV1 = interp1(I1, V1, pI1);
    pV2 = interp1(I2, V2, pI2);

    % calculate Rs
    Rs = (pV2 - pV1) ./ (pI1 - pI2);

    % calculate estimated error
    if Rsh1>0 && Rsh2>0 && nNsVth>0 && Io>0
        delta = log(Rsh1./Rsh2) + log((delI.*Rsh2./nNsVth + log(Io.*Rsh2./nNsVth)) ./ ...
                    (delI.*Rsh1./nNsVth + log(Io.*Rsh1./nNsVth)));
        delta2 = 1.04*( log(delI.*Rsh1 + nNsVth.*log(Io.*Rsh1./nNsVth)) ./ ...
                         (delI.*Rsh1 + nNsVth.*log(Io.*Rsh1./nNsVth)) ...
                      - log(delI.*Rsh2 + nNsVth.*log(Io.*Rsh2./nNsVth)) ./ ...
                           (delI.*Rsh2 + nNsVth.*log(Io.*Rsh2./nNsVth)) );

        errest = nNsVth./(pI1 - pI2).*(delta + delta2);
    else
        errest = NaN;
    end

else
    
    Rs = NaN; errest = NaN;

end


