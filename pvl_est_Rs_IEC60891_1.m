function Rs = pvl_est_Rs_IEC60891_1(IVCurves, aIsc, bVoc)
% PVL_EST_RS_IEC60891_1 estimates Rs using method 1 in IEC 60891.
%
% Syntax
%   Rs = pvl_est_Rs_IEC60891_1(IVCurves, aIsc, bVoc)
%
% Description
%   Method 1 in IEC60891 translates a point on an IV curve measured at 
%   irradiance G1 and cell temperature Tc1 to the corresponding point on an 
%   unobserved IV curve at irradiance G2 and cell temperature Tc2. The 
%   translation reduces the voltage value by a voltage drop across the 
%   ‘internal series resistance’ of the module ([IEC 60891], Eq. 2). A value
%   for the resistance Rs is found by an iterative search over a set of IV 
%   curves measured at different irradiance levels and constant cell 
%   temperature: first, current for each IV curve is translated linearly 
%   to a common irradiance, then Rs is found by minimizing the variance of
%   Pmp of the translated IV curves.
%
%   pvl_est_Rs_IEC60891_1 assumes that the IV curves have been translated 
%   to a common cell temperature.
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
%   Rs - the series resistance value in ohms.
%  
% References
%   [1] IEC60891 Ed. 2 2009. Procedures for temperature and irradiance 
%   corrections to measured I-V characteristics of crystalline silicon 
%   photovoltaic (PV) devices.
%
%   [2] C. Hansen and B. King, "Determining series resistance for
%   equivalent circuit models of a PV module", in 45th IEEE Photovoltaic
%   Specialist Conference, Waikoloa, HI, 2018.

% first check if cell temperatures are all equal
Tc = unique(round([IVCurves(:).Tc], 2));

if length(Tc)>1
    warning('est_Rs_IEC60891_1: IV curves must have equal cell temperature')
    Rs = NaN;
    return
else
    % find maximum irradiance level
    G = max([IVCurves(:).Ee]);
    objfun = @(Rs) max_Pmp_diff(Rs, IVCurves, G, Tc, aIsc, bVoc);
    options = optimset('TolX', 1e-6);
    Rs = fminsearch(objfun, 0.0, options);
    ['Estimated Rs : ' num2str(Rs) ', variance ' num2str(max_Pmp_diff(Rs, IVCurves, G, Tc, aIsc, bVoc))];
end


function val = max_Pmp_diff(Rs, IVCurves, G, Tc, aIsc, bVoc)

% calculates maximum absolute difference in Pmp after translating IVCurves
% to the highest irradiance

for i=1:length(IVCurves)
    tmpIVCurves(i) = pvl_translate_IV_curve_IEC60891_1(IVCurves(i), G, Tc, aIsc, bVoc, Rs, 0);
end
Pmp = [tmpIVCurves(:).Pmp];
val = max(abs((Pmp - median(Pmp))./median(Pmp)));


