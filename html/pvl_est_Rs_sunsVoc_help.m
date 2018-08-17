%% pvl_est_Rs_sunsVoc
% PVL_EST_RS_SUNSVOC Estimate Rs using the suns-Voc method.
%
%% Syntax
%%
% *  |[Rs, errest] = pvl_est_Rs_sunsVoc(sunsV, sunsI, Isc, Vmp, Imp)|
% *  |[Rs, errest] = pvl_est_Rs_sunsVoc(sunsV, sunsI, Isc, Vmp, Imp, Rsh1, Rsh2, nNsVth, Io, Isc2)|
%
%% Description:
%%
% The suns-Voc method estimates series resistance using the suns-Voc 
% curve and a measured IV curve. It is assumed that both the suns-Voc and
% measured IV curves are at cell temperature of 25C. The suns-Voc curve
% is expressed as a pair of vectors (sunsV, sunsI) where sunsV is a
% vector of Voc at irradiance levels E, and sunsI = Isc0 - Isc(E) where
% Isc0 is the short circuit current at STC (1000 W/m2 and 25C). A point
% is selected on the (sunsV, sunsI) curve where the pseudo-current sunsI
% is equal to the measured Imp.
%
% If optional arguments are provided, the difference between the returned
% Rs value and the Rs parameter for the single diode equation is estimated,
% see [2]. For this estimate, shunt resistance and short circuit current
% are required for the IV curve at irradiance E = (1 - Imp / Isc0)*1000 W/m2.
%
%% Inputs:
% * *|sunsV|* - a vector of voltage for suns-Voc curve.
% * *|sunsI|* - a vector of pseudo-current in amps for the suns-Voc curve.
% * *|Isc|* - short circuit current in amps for the measured IV curve.
% * *|Vmp|* - voltage at the maximum power point for the measured IV curve.
% * *|Imp|* - current at the maximum power point for the measured IV curve.
% * *|Rsh1|* - (optional) shunt resistance, in ohms, for the measured IV curve.
% * *|Rsh2|* - (optional) shunt resistance, in ohms, for the IV curve at irradiance
%     E = (1 - Imp / Isc0) * 1000 W/m2 where Isc0 is short circuit current
%     in amps at STC.
% * *|nNsVth|* - (optional) the product n (diode factor) x Ns (cells in series)
%     x Vth (thermal voltage per cell) for both IV curves.
% * *|Io|* - (optional) the dark current, in A, for both IV curves.
% * *|Isc2|* - (optional) short circuit current, in amps, for the IV curve 
%     at irradiance E = (1 - Imp / Isc0) * 1000 W/m2 where Isc0 is short 
%     circuit current in amps at STC.
% 
%% Outputs:
% * *|Rs|* - the series resistance value in ohms.
% * *|errest|* - the estimated difference between Rs and the series resistance
%     parameter for the single diode equation.
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


% Calculate IV curve at STC
Ee = 1000;
Tc = 25;
nPts = 100;
[IL, I0, Rs, Rshf, nNsVth] = pvl_calcparams_desoto(Ee, Tc, param.aIsc, param, EgRef, dEgdT);
IVcurve_STC = pvl_singlediode(IL, I0, Rs, Rshf, nNsVth, nPts);


% simulate suns-Voc curve
suns = 10:10:1000; suns = suns(:);
sunsTc = ones(size(suns)).*25;
sunsIsc0 = IVcurve_STC.Isc;
sunsIsc = suns./1000.*sunsIsc0;

[IL, I0, Rs, Rsh, nNsVth] = pvl_calcparams_desoto(suns, sunsTc, param.aIsc, param, EgRef, dEgdT);
res = pvl_singlediode(IL, I0, Rs, Rsh, nNsVth, nPts);
sunsVoc.V = [res(:).Voc];
sunsVoc.I = (1 - suns(:)/1000)*sunsIsc0;

% Estimate Rs
Rs = pvl_est_Rs_sunsVoc(sunsVoc.V, sunsVoc.I, IVcurve_STC.Isc, IVcurve_STC.Vmp, IVcurve_STC.Imp)


%% References:
% * [1] D. Pysch, A. Mette, S. W. Glunz, “A review and comparison of 
%     different methods to determine the series resistance of solar cells",
%     Solar. Energy Materials and Cells 91, pp. 1698-1706, 2007.
%
% * [2] C. Hansen and B. King, "Determining series resistance for
%     equivalent circuit models of a PV module", in 45th IEEE Photovoltaic
%     Specialist Conference, Waikoloa, HI, 2018.
%
%% See also 
% <pvl_est_Rs_Bowden_help.html |pvl_est_Rs_Bowden|> ,        
% <pvl_est_Rs_Swanson_help.html |pvl_est_Rs_Swanson|> ,        
% <pvl_est_Rs_Pysch_help.html |pvl_est_Rs_Pysch|> ,        
% <pvl_est_Rs_IEC60891_1_help.html |pvl_est_Rs_IEC60891_1|> ,        
% <pvl_est_Rs_IEC60891_2_help.html |pvl_est_Rs_IEC60891_2|>
%
%%
% Copyright 2018 Sandia National Laboratories

