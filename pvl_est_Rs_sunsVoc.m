function [Rs, errest] = pvl_est_Rs_sunsVoc(sunsV, sunsI, Isc, Vmp, Imp, varargin)
% PVL_EST_RS_SUNSVOC estimates Rs using the suns-Voc method.
%
% Description
%   The suns-Voc method estimates series resistance using the suns-Voc 
%   curve and a measured IV curve. It is assumed that both the suns-Voc and
%   measured IV curves are at cell temperature of 25C. The suns-Voc curve
%   is expressed as a pair of vectors (sunsV, sunsI) where sunsV is a
%   vector of Voc at irradiance levels E, and sunsI = Isc0 - Isc(E) where
%   Isc0 is the short circuit current at STC (1000 W/m2 and 25C). A point
%   is selected on the (sunsV, sunsI) curve where the pseudo-current sunsI
%   is equal to the measured Imp.
%
%   If optional arguments are provided, the difference between the returned
%   Rs value and the Rs parameter for the single diode equation is estimated,
%   see [2]. For this estimate, shunt resistance and short circuit current
%   are required for the IV curve at irradiance E = (1 - Imp / Isc0)*1000
%   W/m2.
%
% Syntax
%   [Rs, errest] = pvl_est_Rs_sunsVoc(sunsV, sunsI, Isc, Vmp, Imp)
%   [Rs, errest] = pvl_est_Rs_sunsVoc(sunsV, sunsI, Isc, Vmp, Imp, Rsh1, Rsh2, nNsVth, Io, Isc2)
%
% Inputs:
%   sunsV - a vector of voltage for the suns-Voc curve.
%   sunsI - a vector of pseudo-current in amps for the suns-Voc curve.
%   Isc - short circuit current in amps for the measured IV curve.
%   Vmp - voltage at the maximum power point for the measured IV curve.
%   Imp - current at the maximum power point for the measured IV curve.
%   Rsh1 - (optional) shunt resistance in ohms for the measured IV curve.
%   Rsh2 - (optional) shunt resistance in ohms for the IV curve at irradiance
%     E = (1 - Imp / Isc0) * 1000 W/m2 where Isc0 is short circuit current
%     in amps at STC.
%   nNsVth - (optional) the product n (diode factor) x Ns (cells in series)
%     x Vth (thermal voltage per cell) for both IV curves.
%   Io - (optional) the dark current, in A, for both IV curves.
%   Isc2 - (optional) short circuit current in amps for the IV curve at irradiance
%     E = (1 - Imp / Isc0) * 1000 W/m2 where Isc0 is short circuit current
%     in amps at STC.
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
addRequired(p, 'sunsV', @isvector);
addRequired(p, 'sunsI', @isvector);
addRequired(p, 'Isc', @isnumeric);
addRequired(p, 'Vmp', @isnumeric);
addRequired(p, 'Imp', @isnumeric);

addOptional(p, 'Rsh1', NaN, @isnumeric);
addOptional(p, 'Rsh2', NaN, @isnumeric);
addOptional(p, 'nNsVth', NaN, @isnumeric);
addOptional(p, 'Io', NaN, @isnumeric);
addOptional(p, 'Isc2', NaN, @isnumeric);

parse(p, sunsV, sunsI, Isc, Vmp, Imp, varargin{:});
sunsV = p.Results.sunsV;
sunsI = p.Results.sunsI;
Isc = p.Results.Isc;
Vmp = p.Results.Vmp;
Imp = p.Results.Imp;

Rsh1 = p.Results.Rsh1;
Rsh2 = p.Results.Rsh2;
nNsVth = p.Results.nNsVth;
Io = p.Results.Io;
Isc2 = p.Results.Isc2;


% find voltage on suns-Voc at Imp
pVoc = interp1(sunsI, sunsV, Imp);

if pVoc>0 && ~isnan(pVoc)
    
    Rs = (pVoc - Vmp) ./ Imp;

    if Rsh1>0 && Rsh2>0 && Io>0 && nNsVth>0
        % calculate estimated error
        delta = log(Rsh1./Rsh2) + log((nNsVth.*log(Io.*Rsh2./nNsVth) + Isc2.*Rsh2) ./ ...
                                     (nNsVth.*log(Io.*Rsh1./nNsVth) + (Isc - Imp).*Rsh1));
        delta2 = 1.04 * ( log((nNsVth.*log(Io.*Rsh1./nNsVth) + (Isc - Imp).*Rsh1)) ./ ...
                              (nNsVth.*log(Io.*Rsh1./nNsVth) + (Isc - Imp).*Rsh1) - ...
                          log((nNsVth.*log(Io.*Rsh2./nNsVth) + Isc2.*Rsh2)) ./ ...
                              (nNsVth.*log(Io.*Rsh2./nNsVth) + Isc2.*Rsh2) );

        errest = nNsVth./Imp.*(delta + delta2);

    else
        errest = NaN;
    end
    
else
    
    Rs = NaN; errest = NaN;
    
end



