%% pvl_est_Rs_Pysch
% PVL_EST_RS_PYSCH Estimate Rs using the method of Pysch.
%
%% Syntax
%%
% *  |[Rs, errest] = pvl_est_Rs_Pysch(IVCurves, delI)|
% *  |[Rs, errest] = pvl_est_Rs_Pysch(IVCurves, delI, Rsh, nNsVth, Io)|
%
%% Description:
%%
% The method of Pysch [1] extends the Swanson method to use multiple 
% IV curves. IV curves are assumed to be at different irradiance levels
% but the same cell temperature. A point is selected on each IV curve where 
% the current is less than Isc by the input delI. If optional arguments are
% provided, the difference between the returned Rs value and the Rs
% parameter for the single diode equation is estimated, see [2].
%
%% Inputs:
% * *|IVCurves|* - a structure array for IV curves including fields V and
%     I.
% * *|delI|* - offset from short circuit at which to choose
%     IV curve points.
% * *|Rsh|* - (optional) a vector of shunt resistance, in ohms, for each IV curve.
% * *|nNsVth|* - (optional) a vector containing the product n (diode factor)
%      x Ns (cells in series) x Vth (thermal voltage per cell) for each IV curve.
% * *|Io|* - (optional) a vector of the dark current, in A, for both IV curves.
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

% Calculate set of IV curves
Ee = 900:20:1100;
Tc = 25;
nPts = 100;
for i=1:length(Ee)
    [IL, I0, Rs, Rsh1, nNsVth] = pvl_calcparams_desoto(Ee(i), Tc, param.aIsc, param, EgRef, dEgdT);
    IVcurves(i) = pvl_singlediode(IL, I0, Rs, Rsh1, nNsVth, nPts);
end

% Current differential
STCidx = find(Ee==1000);

delI = 0.95*(IVcurves(STCidx).Isc - IVcurves(STCidx).Imp);

% Estimate Rs
Rs = pvl_est_Rs_Pysch(IVcurves, delI)



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
% <pvl_est_Rs_sunsVoc_help.html |pvl_est_Rs_sunsVoc|> ,        
% <pvl_est_Rs_IEC60891_1_help.html |pvl_est_Rs_IEC60891_1|> ,        
% <pvl_est_Rs_IEC60891_2_help.html |pvl_est_Rs_IEC60891_2|>
%
%%
% Copyright 2018 Sandia National Laboratories

