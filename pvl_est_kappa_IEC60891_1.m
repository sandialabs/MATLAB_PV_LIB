function kappa = pvl_est_kappa_IEC60891_1(IVCurves, aIsc, bVoc)
% PVL_EST_KAPPA_IEC60891_1 estimates the curve correction factor kappa for 
% IEC 60891 method 1.
%
% Syntax
%   kappa = pvl_est_kappa_IEC60891_1(IVCurves, aIsc, bVoc)
%
% Description
%   Methods in IEC60891 translates a point on an IV curve measured at 
%   irradiance G1 and cell temperature Tc1 to the corresponding point on an 
%   unobserved IV curve at irradiance G2 and cell temperature Tc2. The 
%   translation reduces the voltage value by a voltage drop across the 
%   ‘internal series resistance’ of the module ([IEC 60891], Eq. 2). A value
%   for the curve correction factor interpreted as the temperature coefficient
%   of the series resistance is found by interative search using a set of 
%   IV curves measured at different cell temperatures and at a common irradiance:
%   first, voltage and current for each IV curve are translated to a common
%   temperature, then kappa is found by minimizing the variance of
%   Pmp of the translated IV curves.
%
%   pvl_est_kappa_IEC60891_1 assumes that the IV curves have been translated 
%   to a common irradiance.
%
% Input:
%   IVCurves - A structure array with the following fields:
%      IVCurves.Isc - short circuit current in amperes.
%      IVCurves.Voc - open circuit voltage in volts.
%      IVCurves.Imp - current at maximum power point in amperes. 
%      IVCurves.Vmp - voltage at maximum power point in volts.
%      IVCurves.Pmp - power at maximum power point in watts.
%      IVCurves.V - vector of voltage in volts. 
%      IVCurves.I - vector of current in amperes.
%      IVCurves.Ee - Effective irradiance (W/m2)
%      IVCurves.Tc - cell temperature (C)
%   aIsc - temperature coefficient for short circuit current in A/C
%   bVoc - temperature coefficient for open circuit voltage in V/C
%
% Output:
%   kappa - the curve correction factor in ohm/K.
%  
% References
%   [1] IEC60891 Ed. 2 2009. Procedures for temperature and irradiance 
%   corrections to measured I-V characteristics of crystalline silicon 
%   photovoltaic (PV) devices.
%
%   [2] C. Hansen and B. King, "Determining series resistance for
%   equivalent circuit models of a PV module", in 45th IEEE Photovoltaic
%   Specialist Conference, Waikoloa, HI, 2018.

% first check if irradiance levels are all equal
Ee = unique(round([IVCurves(:).Ee], 2));

if length(Ee)>1
    warning('est_kappa_IEC60891_1: IV curves must have equal irradiance')
    kappa = NaN;
    return
else
    % find minimum cell temperature
    Tc_min = max([IVCurves(:).Tc]);
    objfun = @(kappa) max_Pmp_diff(kappa, IVCurves, Ee, Tc_min, aIsc, bVoc);
    options = optimset('TolX', 1e-6);
    kappa = fminsearch(objfun, 0.0, options);
    ['Estimated kappa : ' num2str(kappa) ', variance ' num2str(max_Pmp_diff(kappa, IVCurves, Ee, Tc_min, aIsc, bVoc))];
end


function val = max_Pmp_diff(kappa, IVCurves, G, Tc, aIsc, bVoc)

% calculates maximum absolute difference in Pmp after translating IVCurves
% to the highest irradiance

for i=1:length(IVCurves)
    tmpIVCurves(i) = pvl_translate_IV_curve_IEC60891_1(IVCurves(i), G, Tc, aIsc, bVoc, 0, kappa);
end
Pmp = [tmpIVCurves(:).Pmp];
val = max(abs((Pmp - median(Pmp))./median(Pmp)));


