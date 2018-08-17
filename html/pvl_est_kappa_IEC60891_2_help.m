%% pvl_est_kappa_IEC60891_2
% PVL_EST_KAPPA_IEC60891_2 Estimate the curve correction factor kappa for 
% IEC 60891 method 2.
%
%% Syntax
%%
% |kappa = pvl_est_kappa_IEC60891_2(IVCurves, aIsc, bVoc, Rs)|
%
%% Description:
%%
% Methods in IEC60891 translates a point on an IV curve measured at 
% irradiance G1 and cell temperature Tc1 to the corresponding point on an 
% unobserved IV curve at irradiance G2 and cell temperature Tc2. The 
% translation reduces the voltage value by a voltage drop across the 
% ‘internal series resistance’ of the module ([IEC 60891], Eq. 2). A value
% for the curve correction factor interpreted as the temperature coefficient
% of the series resistance is found by interative search using a set of 
% IV curves measured at different cell temperatures and at a common irradiance:
% first, voltage and current for each IV curve are translated to a common
% temperature, then kappa is found by minimizing the variance of
% Pmp of the translated IV curves.
%
% |pvl_est_kappa_IEC60891_2| assumes that the IV curves have been translated 
% to a common irradiance.
%
%% Inputs:
% * *|IVCurves|* - A structure array with the following fields: 
% * * *|IVCurves.Isc|* - short circuit current in amperes.
% * * *|IVCurves.Voc|* - open circuit voltage in volts.
% * * *|IVCurves.Imp|* - current at maximum power point in amperes. 
% * * *|IVCurves.Vmp|* - voltage at maximum power point in volts.
% * * *|IVCurves.Pmp|* - power at maximum power point in watts.
% * * *|IVCurves.V|* - vector of voltage in volts. 
% * * *|IVCurves.I|* - vector of current in amperes.
% * * *|IVCurves.Ee|* - Effective irradiance (W/m2).
% * * *|IVCurves.Tc|* - cell temperature (C).
% * *|aIsc|* - relative temperature coefficient for short circuit current 
%      in 1/C (not A/C). 
% * *|bVoc|* - relative temperature coefficient for open circuit voltage 
%      in 1/C (not V/C).
% * *|Rs|* - the series resistance estimated for curve translation in ohms, see
%     |pvl_est_Rs_IEC60891_2|
% 
%% Outputs:
% * *|kappa|* - the curve correction factor in ohm/K.
%
%% Example
clearvars

% Set up parameters for a representative 60 cell cSi module using the Desoto model
Const.q = 1.60218E-19;
Const.k = 1.38066E-23;
Const.E0 = 1000;
Const.T0 = 25;

param.aIsc = 0.0008;  % A/C
param.bVoc = -0.1900; % V/C
 
param.Rs_ref = 0.2;
param.Rsh_ref = 1000;
param.IL_ref = 8.0;
param.I0_ref = 5e-10;
 
param.a_ref = 1.05 * 60 * Const.k/Const.q * (273.15 + Const.T0);
 
EgRef = 1.121;
dEgdT = -0.0002677;

% Calculate set of IV curves
Ee = 1000;
Tc = 25:1:35;
nPts = 100;
for i=1:length(Tc)
    [IL, I0, Rs, Rsh1, nNsVth] = pvl_calcparams_desoto(Ee, Tc(i), param.aIsc, param, EgRef, dEgdT);
    IVcurves(i) = pvl_singlediode(IL, I0, Rs, Rsh1, nNsVth, nPts);
end

% Add Ee and Tc to IVcurves
for i=1:length(Tc)
    IVcurves(i).Ee = Ee;
    IVcurves(i).Tc = Tc(i);
end

% Calculate relative temperature coefficients

STC_Isc = IVcurves(6).Isc;
STC_Voc = IVcurves(6).Voc;
aIsc_rel = param.aIsc / STC_Isc;
bVoc_rel = param.bVoc / STC_Voc;

% Use known Rs for the curve correction factor. Otherwise, estimate Rs from
% a set of IV curves with fixed temperature and varying irradiance, using
% pvl_est_Rs_IEC60891_2
Rs_est = param.Rs_ref;

% Estimate kappa
kappa = pvl_est_kappa_IEC60891_2(IVcurves, aIsc_rel, bVoc_rel, Rs_est)


%% References:
% * [1] IEC60891 Ed. 2 2009. Procedures for temperature and irradiance 
%     corrections to measured I-V characteristics of crystalline silicon 
%     photovoltaic (PV) devices.
%
% * [2] C. Hansen and B. King, "Determining series resistance for
%     equivalent circuit models of a PV module", in 45th IEEE Photovoltaic
%     Specialist Conference, Waikoloa, HI, 2018.
%
%% See also 
% <pvl_est_Rs_IEC60891_2_help.html |pvl_est_Rs_IEC60891_2|>,
% <pvl_translate_IV_curve_IEC60891_2_help.html |pvl_translate_IV_curve_IEC60891_2|>
%
%%
% Copyright 2018 Sandia National Laboratories

