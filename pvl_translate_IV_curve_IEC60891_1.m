function Result = pvl_translate_IV_curve_IEC60891_1(IVCurves, G, Tc, aIsc, bVoc, Rs, kappa)
% PVL_TRANSLATE_IV_CURVE_IEC60891_1 translates IV curves to a target 
% irradiance G and temperature Tc using method 1 in IEC 60891.
%
% Syntax
%   Result = pvl_translate_IV_curve_IEC60891_1(IVCurves, G, Tc, aIsc, bVoc, Rs, kappa)
%
% Input:
%   IVCurves - A structure with the following fields:
%      IVCurves.Isc - Column vector of short circuit current in amperes.
%      IVCurves.Voc - Column vector of open circuit voltage in volts.
%      IVCurves.Imp - Column vector of current at maximum power point in amperes. 
%      IVCurves.Vmp - Column vector of voltage at maximum power point in volts.
%      IVCurves.Pmp - Column vector of power at maximum power point in watts.
%      IVCurves.V - Array of voltages in volts. Row n corresponds to IV 
%        curve n, with V=0 in the leftmost column and V=Voc in the rightmost
%        column.
%      IVCurves.I - Array of currents in amperes. Row n corresponds to IV 
%        curve n, with I=Isc in the leftmost column and I=0 in the 
%        rightmost column.
%      IVCurves.Ee - column vector of effective irradiance (W/m2)
%      IVCurves.Tc - column vector of cell temperature (C)
%   G - the target irradiance (W/m2)
%   Tc - the target cell temperature (C)
%   aIsc - temperature coefficient for short circuit current in A/C
%   bVoc - temperature coefficient for open circuit voltage in V/C
%   Rs - internal series resistance parameter for the translation, in ohms.
%   kappa - curve correction factor interpreted as the temperature
%     coefficient of the internal series resistance, in ohms/C.
%
% Output:
%   Result - a structure with the same fields as IVCurves, containing the
%   translated IV curves
%
% References
%   [1] IEC60891 Ed. 2 2009. Procedures for temperature and irradiance 
%   corrections to measured I-V characteristics of crystalline silicon 
%   photovoltaic (PV) devices.

N = size(IVCurves.V, 1);

clearvars Result;

for i=1:N
    I = IVCurves.I(i, :);
    V = IVCurves.V(i, :);
    
    tI = I + IVCurves.Isc(i)*(G/IVCurves.Ee(i) - 1) + aIsc*(Tc - IVCurves.Tc(i));
    tV = V - Rs*(tI - I) - kappa*tI*(Tc - IVCurves.Tc(i)) + bVoc*(Tc - IVCurves.Tc(i));
    
    Result.Ee(i) = G;
    Result.Tc(i) = Tc;
    tP = tI.*tV;
    [tPmp, idx] = max(tP);
    Result.Pmp(i) = tPmp;
    Result.Vmp(i) = tV(idx);
    Result.Imp(i) = tI(idx);
    Result.Isc(i) = interp1(tV, tI, 0, 'linear', 'extrap');
    Result.Voc(i) = interp1(tI, tV, 0, 'linear', 'extrap');
    Result.I(i,:) = tI;
    Result.V(i,:) = tV;
    
end


    