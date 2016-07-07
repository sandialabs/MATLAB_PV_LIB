%% pvl_PVsyst_parameter_estimation
% Estimate parameters for the PVsyst v6 module performance model.
%
%% Syntax
% * |[PVsyst oflag] = pvl_PVsyst_parameter_estimation(IVCurves, Specs,
% Const, maxiter, eps1, graphic)|
%
%% Description
% |pvl_PVsyst_parameter_estimation| estimates parameters for the PVsyst 
% version 6 module performance model [1,2,3]. Estimation methods are
% documented in [4, 5, 6].
%
%% Inputs
% * *|IVCurves|* - a structure array ontaining IV curve data in the following
% fields:
% * *|IVCurves(i).I|* - vector of current (A) (same length as V)
% * *|IVCurves(i).V|* - vector of voltage (V) (same length as I)
% * *|IVCurves(i).Ee|* - effective irradiance (W/m^2), i.e., POA broadband
%       irradiance adjusted by spectral mismatch modifier
% * *|IVCurves(i).Tc|* - cell temperature (C)
% * *|IVCurves(i).Isc|* - short-circut current of IV curve (A)
% * *|IVCurves(i).Voc|* - open-curcut voltage of IV curve (V)
% * *|IVCurves(i).Imp|* - current at max power point of IV curve (A)
% * *|IVCurves(i).Vmp|* - voltage at max power point of IV curve (V)
%
% * *|Specs|* - a structure containing module-level values
% * *|Specs.Ns|* - number of cells in series
% * *|Specs.aIsc|* - temperature coefficient of Isc (A/C)
%
% * *|Const|* - a structure containing physical and other constants
% * *|Const.E0|* - effective irradiance at STC, normally 1000 W/m2
% * *|Const.T0|* - cell temperature at STC, normally 25 C
% * *|Const.k|* - 1.38066E-23 J/K (Boltzmann's constant)
% * *|Const.q|* - 1.60218E-19 Coulomb (elementary charge)
%
%% Optional inputs
% * *|maxiter|* - an integer setting the maximum number of iterations for the 
%     parameter updating part of the algorithm. Default value is 5
%
% * *|eps1|* - the desired tolerance for convergence for the IV curve fitting.
%     The iterative parameter updating stops when absolute values of the
%     relative change in mean, max and standard deviation of Imp, Vmp and Pmp
%     between iterations are all less than eps1, or when the number of 
%     iterations exceeds maxiter.  Default value of eps1 is 1e-3 (0.0001%). 
%
% * *|graphic|* - a boolean, if true then plots are produced during the 
%     parameter estimation process. Default is false
%
%% Outputs
% * *|PVsyst|* - a structure containing the model parameters
% * *|PVsyst.IL_ref|* - light current (A) at STC
% * *|PVsyst.Io_ref|* - dark current (A) at STC
% * *|PVsyst.eG|* - effective band gap (eV) at STC
% * *|PVsyst.Rsh_ref|* - shunt resistance (ohms) at STC
% * *|PVsyst.Rsh0|* - shunt resistance (ohms) at zero irradiance 
% * *|PVsyst.Rshexp|* - exponential factor defining decrease in 
%      Rsh with increasing effective irradiance
% * *|PVsyst.Rs_ref|* - series resistance (ohms) at STC
% * *|PVsyst.gamma_ref|* - diode (ideality) factor at STC
% * *|PVsyst.mugamma|* - temperature coefficient for diode (ideality) factor
% * *|PVsyst.Iph|* - vector of values of light current Iph estimated for each IV
%    curve
% * *|PVsyst.Io|* - vector of values of dark current Io estimated for each IV
%       curve
% * *|PVsyst.Rsh|* - vector of values of shunt resistance Rsh estimated for each IV
%    curve
% * *|PVsyst.Rs|* - vector of values of series resistance Rs estimated for each IV
%       curve
% * *|PVsyst.u|* - logical array indicating IV curves with parameter values deemed 
%       reasonable by the private function filter_params
%
% * *|oflag|* - Boolean indicating success or failure of estimation of the
%    diode (ideality) factor parameter.  If failure, then no parameter values
%    are returned.
%
%% Example

clearvars

% load IV curve data for a 36 cell Mitsubishi cSi module
load 'PVsyst_demo.mat'

% Build structure for constants
Const.E0 = 1000; % W/m2
Const.T0 = 25; % C
Const.k = 1.38066E-23; % J/K
Const.q = 1.60218E-19; % c


% set control variables for parameter algorithm
maxiter = 20;
eps1 = NaN; % use default
graphic = false;

% Estimate PVsyst model parameters. Structure PVsyst contains the model
% parameters.  Logical oflag indicates success or failure. Logical array u
% is true for IV curves where valid parameter sets were determined.
% Remaining output variables contain parameter values for each IV curve.
[PVsyst oflag] = pvl_PVsyst_parameter_estimation(IVCurves, Specs, Const, maxiter, eps1, graphic);

% Calculate 5 parameters for each IV curve using PVsyst model
[IL, Io, Rs, Rsh, nNsVth] = pvl_calcparams_PVsyst([IVCurves.Ee],[IVCurves.Tc],Specs.aIsc,PVsyst);
% Compute IV curve points using 5 parameters
Modeled = pvl_singlediode(IL, Io, Rs, Rsh, nNsVth);

% Compute IV curve at STC
[IL0, Io0, Rs0, Rsh0, nNsVth0] = pvl_calcparams_PVsyst(1000,25,Specs.aIsc,PVsyst);
% Compute IV curve points using 5 parameters. Structure STCModeled includes
% fields Imp, Vmp, Isc, Voc, Pmp which will be modeled values at STC
STCModeled = pvl_singlediode(IL0, Io0, Rs0, Rsh0, nNsVth0);

FF0 = (STCModeled.Imp*STCModeled.Vmp)/(STCModeled.Isc*STCModeled.Voc);  % Fill factor at STC

% Calculate temperature coefficients for Vmp, Imp and Voc from simulated IV
% curves at Ee = 1000. Create vector of cell temperatures for simulation to determine temperature
% coefficients
vTc = [15 25 35 45 55 65]';
% Compute sets of parameter values for simulation to determine temperature
% coefficients
[ILtc, Iotc, Rstc, Rshtc, nNsVthtc] = pvl_calcparams_PVsyst(1000,vTc,Specs.aIsc,PVsyst);
% Compute IV curves for simulation to determine temperature
% coefficients
Modeledtc = pvl_singlediode(ILtc, Iotc, Rstc, Rshtc, nNsVthtc);

% extract IV curve points and estimate temperature coefficients by
% regression
col1 = ones(size(vTc));
[beta_Vmp]=[vTc-25 col1]\[Modeledtc.Vmp];
betaVmp = beta_Vmp(1);      % temperature coefficient (V/C) for voltage at maximum power
[beta_Voc]=[vTc-25 col1]\[Modeledtc.Voc];
betaVoc = beta_Voc(1);      % temperature coefficient (V/C) for open circuit voltage
[alpha_Imp]=[vTc-25 col1]\[Modeledtc.Imp];
alphaImp = alpha_Imp(1);    % temperature coefficient (W/C) for current at maximum power

gammaPmp = alphaImp*STCModeled.Vmp+betaVmp*STCModeled.Imp;  % temperature coefficient for maximum power


% Compare calculated and measured IV curve points

Data.Isc = [IVCurves.Isc];
Data.Imp = [IVCurves.Imp];
Data.Voc = [IVCurves.Voc];
Data.Vmp = [IVCurves.Vmp];
Data.Ee = [IVCurves.Ee];
Data.Isc = Data.Isc(:);
Data.Imp = Data.Imp(:);
Data.Voc = Data.Voc(:);
Data.Vmp = Data.Vmp(:);
Data.Ee = Data.Ee(:);

% Plot comparison
figure('Position',[50 50 50+9*96 50+6*96])

Meas = [Data.Isc, Data.Imp, Data.Imp.*Data.Vmp, Data.Voc, Data.Vmp];
Mod = [Modeled.Isc, Modeled.Imp, Modeled.Pmp, Modeled.Voc, Modeled.Vmp];
Title = {'Isc', 'Imp', 'Pmp', 'Voc', 'Vmp'};

for i = [1:5]
    subplot(2,3,i)
    scatter(Meas(:,i), Mod(:,i), 5, 'k', 'filled')
    hold on
    title(Title(i),'fontsize',12,'fontweight','b')
    xlabel('Measured','fontsize',12)
    ylabel('Modeled','fontsize',12)
    axis square
    datamin = min([Meas(:,i);Mod(:,i)]);
    datamax = max([Meas(:,i);Mod(:,i)]);
    datamin = max(datamin - (datamax-datamin)/10,0);
    datamax = datamax + (datamax-datamin)/10;
    plot([datamin, datamax], [datamin, datamax],'r', 'LineWidth',2)
    axis([datamin datamax datamin  datamax])
    box on

end

%% Sources:
% [1] K. Sauer, T. Roessler, C. W. Hansen, Modeling the Irradiance and 
%     Temperature Dependence of Photovoltaic Modules in PVsyst, 
%     IEEE Journal of Photovoltaics v5(1), January 2015.
%
% [2] A. Mermoud, PV modules modelling, Presentation at the 2nd PV
%     Performance Modeling Workshop, Santa Clara, CA, May 2013
%
% [3] A. Mermoud, T. Lejeune, Performance Assessment of a Simulation Model
%     for PV modules of any available technology, 25th European Photovoltaic
%     Solar Energy Conference, Valencia, Spain, Sept. 2010
%
% [4] C. Hansen, Estimating Parameters for the PVsyst Version 6 Photovoltaic
%     Module Performance Model, Sandia National Laboratories Report
%     SAND2015-8598
%
% [5] C. Hansen, Parameter Estimation for Single Diode Models of 
%     Photovoltaic Modules, Sandia National Laboratories Report
%     SAND2015-2065
%
% [6] C. Hansen, Estimation of Parameters for Single Diode Models using
%     Measured IV Curves, Proc. of the 39th IEEE PVSC, June 2013.


%% See also 
% <pvl_calcparams_PVsyst_help.html |pvl_calcparams_PVsyst|>,  
% <pvl_singlediode_help.html |pvl_singlediode|>  
%%
% Copyright 2015 Sandia National Laboratories



