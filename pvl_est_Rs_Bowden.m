function [Rs, errest] = pvl_est_Rs_Bowden(Vocs, Iscs, V, I, varargin)
% PVL_EST_RS_BOWDEN Estimate Rs using the method of Bowden and Rohatgi [1].
%
% Syntax
% [Rs, errest] = pvl_est_Rs_Bowden(Vocs, Iscs, V, I)
% [Rs, errest] = pvl_est_Rs_Bowden(Vocs, Iscs, V, I, Rshf, Rshs, Io, nNsVth)
%
% Description
%   The method of Bowden and Rohatgi estimates a value for series
%   resistance using values from two IV curves, one under 'shaded'
%   conditions, anticipated to be ~100 W/m2, and another at 'full-sun'
%   conditions, anticipated to be ~1000 W/m2. The cell temperature is 
%   assumed to be the same for both IV curves. Voc and Isc from the shaded 
%   IV curve are used with Isc from the full sun IV curve to locate the
%   desired point on the full sun IV curve.  If optional arguments are
%   provided, the difference between the returned Rs value and the Rs
%   parameter for the single diode equation is estimated, see [2].
%
% Inputs:
%   Vocs - open circuit voltage on the shaded (low irradiance) IV curve, 
%     in V.
%   Iscs - short circuit current on the shaded (low irradiance) IV curve, 
%     in A.
%   V - a vector of voltages for the full sun IV curve. It is assumed that
%     V(1) = 0.
%   I - a vector of currents for the full sun IV curve. It is assumed that
%     I(1) = Isc.
%   Rshs - (optional) shunt resistance, in ohms, for the shaded IV curve.
%   Rshf - (optional) shunt resistance, in ohms, for the full-sun IV curve.
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
%   [1] S. Bowden and Rohatgi, A., “Rapid and Accurate Determination of 
%   Series Resistance and Fill Factor Losses in Industrial Silicon Solar 
%   Cells”, in 17th European Photovoltaic Solar Energy Conference, Munich,
%   Germany, 2001.
%
%   [2] C. Hansen and B. King, "Determining series resistance for
%   equivalent circuit models of a PV module", in 45th IEEE Photovoltaic
%   Specialist Conference, Waikoloa, HI, 2018.

p = inputParser;
addRequired(p, 'Vocs', @isnumeric);
addRequired(p, 'Iscs', @isnumeric);
addRequired(p, 'V', @isvector);
addRequired(p, 'I', @isvector);
addOptional(p, 'Rshs', NaN, @isnumeric);
addOptional(p, 'Rshf', NaN, @isnumeric);
addOptional(p, 'nNsVth', NaN, @isnumeric);
addOptional(p, 'Io', NaN, @isnumeric);

parse(p, Vocs, Iscs, V, I, varargin{:});
Vocs = p.Results.Vocs;
Iscs = p.Results.Iscs;
V = p.Results.V;
I = p.Results.I;
Rshs = p.Results.Rshs;
Rshf = p.Results.Rshf;
nNsVth = p.Results.nNsVth;
Io = p.Results.Io;

% calculate Rs
% find voltage on full-sun IV curve
Iscf = I(1);
Va = interp1(I, V, Iscf - Iscs);

Rs = (Vocs - Va) ./ (Iscf - Iscs);

% estimate difference between Rs and the Rs parameter for the single diode
% equation
if Rshf>0 && Rshs>0 && Io>0 && nNsVth>0
    % calculate estimated error
    deltaest = log(Rshf./Rshs) + log( (Iscs.*Rshs./nNsVth + log(Io.*Rshs./nNsVth)) ...
        ./ (Iscs.*Rshf./nNsVth + log(Io.*Rshf./nNsVth)));
    deltaest2 = 1.04*( log(log(Io.*Rshf./nNsVth) + Iscs.*Rshf./nNsVth) ./ ...
                       (log(Io.*Rshf./nNsVth) + Iscs.*Rshf./nNsVth) ...
                     - log(log(Io.*Rshs./nNsVth) + Iscs.*Rshs./nNsVth) ./ ...
                       (log(Io.*Rshs./nNsVth) + Iscs.*Rshs./nNsVth));

    errest = nNsVth./(Iscf - Iscs).*(deltaest + deltaest2);
else
    errest = NaN;
end



