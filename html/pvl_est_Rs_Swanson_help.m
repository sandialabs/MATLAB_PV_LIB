%% pvl_est_Rs_Swanson
% PVL_EST_RS_SWANSON Estimate Rs using the method of Swanson.
%
%% Syntax
%%
% *  |[Rs, errest] = pvl_est_Rs_Swanson(V1, I1, V2, I2, delI)|
% *  |[Rs, errest] = pvl_est_Rs_Swanson(V1, I1, V2, I2, delI, Rsh1, Rsh2, nNsVth, Io)|
%
%% Description:
%%
% The method of Swanson (see [1]) estimates series resistance using 
% points from two IV curves at unequal but similar irradiance, e.g.,
% at 950 W/m2 and 1000 W/m2. The cell temperature is assumed to be the 
% same for both IV curves. A point is selected on each IV curve where 
% the current is less than Isc by the input delI. If optional arguments are
% provided, the difference between the returned Rs value and the Rs
% parameter for the single diode equation is estimated, see [2].
%
%% Inputs:
% * *|V1|* - a vector of voltage for the first IV curve, in V.
% * *|I1|* - a vector of current for the first IV curve, in A.
% * *|V2|* - a vector of voltage for the second IV curve, in V.
% * *|I2|* - a vector of current for the second IV curve, in A.
% * *|delI|* - offset from short circuit at which to choose
%     IV curve points.
% * *|Rsh1|* - (optional) shunt resistance, in ohms, for the first IV curve.
% * *|Rsh2|* - (optional) shunt resistance, in ohms, for the second IV curve.
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

% Calculate first IV curve
Ee = 1000;
Tc = 25;
nPts = 100;
[IL, I0, Rs, Rsh1, nNsVth] = pvl_calcparams_desoto(Ee, Tc, param.aIsc, param, EgRef, dEgdT);
IVcurve1 = pvl_singlediode(IL, I0, Rs, Rsh1, nNsVth, nPts);

% Calculate second IV curve
Ee = 950;
[IL, I0, Rs, Rsh2, nNsVth] = pvl_calcparams_desoto(Ee, Tc, param.aIsc, param, EgRef, dEgdT);
IVcurve2 = pvl_singlediode(IL, I0, Rs, Rsh2, nNsVth, nPts);

% Current differential
delI = 0.95*(IVcurve1.Isc - IVcurve1.Imp);

% Estimate Rs
Rs = pvl_est_Rs_Swanson(IVcurve1.V, IVcurve1.I, IVcurve2.V, IVcurve2.I, delI)


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
% <pvl_est_Rs_sunsVoc_help.html |pvl_est_Rs_sunsVoc|> ,        
% <pvl_est_Rs_Pysch_help.html |pvl_est_Rs_Pysch|> ,        
% <pvl_est_Rs_IEC60891_1_help.html |pvl_est_Rs_IEC60891_1|> ,        
% <pvl_est_Rs_IEC60891_2_help.html |pvl_est_Rs_IEC60891_2|>
%
%%
% Copyright 2018 Sandia National Laboratories

