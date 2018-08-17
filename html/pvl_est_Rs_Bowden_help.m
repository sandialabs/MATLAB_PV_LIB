%% pvl_est_Rs_Bowden
% PVL_EST_RS_BOWDEN Estimate Rs using the method of Bowden and Rohatgi.
%
%% Syntax
%%
% *  |[Rs, errest] = pvl_est_Rs_Bowden(Vocs, Iscs, V, I)|
% *  |[Rs, errest] = pvl_est_Rs_Bowden(Vocs, Iscs, V, I, Rshf, Rshs, Io, nNsVth)|
%
%% Description:
%%
% The method of Bowden and Rohatgi estimates a value for series
% resistance using values from two IV curves, one under 'shaded'
% conditions, anticipated to be ~100 W/m2, and another at 'full-sun'
% conditions, anticipated to be ~1000 W/m2. The cell temperature is 
% assumed to be the same for both IV curves. Voc and Isc from the shaded 
% IV curve are used with Isc from the full sun IV curve to locate the
% desired point on the full sun IV curve.  If optional arguments are
% provided, the difference between the returned Rs value and the Rs
% parameter for the single diode equation is estimated, see [2].
%
%% Inputs:
% * *|Vocs|* - open circuit voltage on the shaded (low irradiance) IV curve, 
%     in V.
% * *|Iscs|* - short circuit current on the shaded (low irradiance) IV curve, 
%     in A.
% * *|V|* - a vector of voltages for the full sun IV curve. It is assumed that
%     V(1) = 0.
% * *|I|* - a vector of currents for the full sun IV curve. It is assumed that
%     I(1) = Isc.
% * *|Rshs|* - (optional) shunt resistance, in ohms, for the shaded IV curve.
% * *|Rshf|* - (optional) shunt resistance, in ohms, for the full-sun IV curve.
% * *|nNsVth|* - (optional) the product n (diode factor) x Ns (cells in series)
%     x Vth (thermal voltage per cell) for both IV curves.
% * *|Io|* - (optional) the dark current, in A, for both IV curves.
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

% Calculate full sun IV curve
Ee = 1000;
Tc = 25;
nPts = 100;
[IL, I0, Rs, Rshf, nNsVth] = pvl_calcparams_desoto(Ee, Tc, param.aIsc, param, EgRef, dEgdT);
fullsunIVcurve = pvl_singlediode(IL, I0, Rs, Rshf, nNsVth, nPts);

% Calculate shaded IV curve
Ee = 200;
[IL, I0, Rs, Rshs, nNsVth] = pvl_calcparams_desoto(Ee, Tc, param.aIsc, param, EgRef, dEgdT);
shadedIVcurve = pvl_singlediode(IL, I0, Rs, Rshs, nNsVth, nPts);

% Estimate Rs
Rs = pvl_est_Rs_Bowden(shadedIVcurve.Voc, shadedIVcurve.Isc, fullsunIVcurve.V, fullsunIVcurve.I)


%% References:
% * [1] S. Bowden and Rohatgi, A., “Rapid and Accurate Determination of 
%     Series Resistance and Fill Factor Losses in Industrial Silicon Solar 
%     Cells”, in 17th European Photovoltaic Solar Energy Conference, Munich,
%     Germany, 2001.
%
% * [2] C. Hansen and B. King, "Determining series resistance for
%     equivalent circuit models of a PV module", in 45th IEEE Photovoltaic
%     Specialist Conference, Waikoloa, HI, 2018.
%
%% See also 
% <pvl_est_Rs_Swanson_help.html |pvl_est_Rs_Swanson|> ,        
% <pvl_est_Rs_sunsVoc_help.html |pvl_est_Rs_sunsVoc|> ,        
% <pvl_est_Rs_Pysch_help.html |pvl_est_Rs_Pysch|> ,        
% <pvl_est_Rs_IEC60891_1_help.html |pvl_est_Rs_IEC60891_1|> ,        
% <pvl_est_Rs_IEC60891_2_help.html |pvl_est_Rs_IEC60891_2|>
%
%%
% Copyright 2018 Sandia National Laboratories

