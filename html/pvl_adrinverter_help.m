%% pvl_adrinverter
% Anton Driesse's Grid-Connected PV Inverter Model
%
%% Syntax
% * |ACPower = pvl_adrinverter(Inverter, Vdc, Pdc)| 
%% Description
% Computes AC power output using Anton Driesse's Grid-Connected PV Inverter model.
%
%% Inputs
%%
% * *|Vdc|* - a vector of modeled DC voltage
% * *|Pdc|* - a vector of modeled DC power
% * *|Inverter|* - A struct defining the inverter to be used, giving the
%     inverter performance parameters according to the Sandia
%     Grid-Connected Photovoltaic Inverter Model (SAND 2007-5036) [1]. A set of
%     inverter performance parameters are provided with PV_LIB (|\Required Data\SandiaInverterDatabaseSAM2012.11.30.mat|),
%     or may be generated from a System Advisor Model (SAM) [2] library using the SAM
%     library reader function <pvl_SAMLibraryReader_SNLInverters_help.html |pvl_SAMLibraryReader_SNLInverters|>. 
%     Required struct components are:
% * *|Inverter.Pac0|* - AC-power output from inverter based on input power and voltage, (W) 
% * *|Inverter.Pdc0|* - DC-power input to inverter, typically assumed to be equal to the PV array maximum
% power, (W)
% * *|Inverter.Vdc0|* - DC-voltage level at which the AC-power rating is achieved at the reference operating
% condition, (V)
% * *|Inverter.Ps0|* - DC-power required to start the inversion process, or self-consumption by inverter,
% strongly influences inverter efficiency at low power levels, (W)
% * *|Inverter.C0|* - Parameter defining the curvature (parabolic) of the relationship between ac-power and
% dc-power at the reference operating condition, default value of zero gives a linear
% relationship, (1/W)
% * *|Inverter.C1|* - Empirical coefficient allowing Pdco to vary linearly with dc-voltage input, default value
% is zero, (1/V)
% * *|Inverter.C2|* - empirical coefficient allowing Pso to vary linearly with dc-voltage input, default value
% is zero, (1/V)
% * *|Inverter.C3|* - empirical coefficient allowing Co to vary linearly with dc-voltage input, default value is
% zero, (1/V)
% * *|Inverter.Pnt|* - ac-power consumed by inverter at night (night tare) to maintain circuitry required to
% sense PV array voltage, (W)
%
%% Outputs
%%
% * *|ACPower|* - a column vector of modeled AC power output given the input 
%     DC voltage, Vdc, and input DC power, Pdc. When ACPower would be 
%     greater than Pmax, it is set to Pmax to represent inverter 
%     "clipping". When ACPower would be less than -Pnt (energy consumed rather
%     than produced) then ACPower is set to -Pnt to represent nightly 
%     power losses. ACPower is not adjusted for maximum power point
%     tracking (MPPT) voltage windows or maximum current limits of the
%     inverter.

%% Example 1
load 'DriesseInverterDatabaseSAM2013.10.mat';
% PV Powered PVP2500 is entry #377
Inverter = ADRInverterDB(377)
%%
% Assume DC power is 1,000 W and voltage is 400 V
%
Pdc = 1000; %DC power is 1000 W
Vdc = 400; % DC violtage is 450 V
ACPower = pvl_adrinverter(Inverter, Vdc, Pdc)
%%
% Inverter efficiency can be calculated as:
Inverter_Efficiency = ACPower/Pdc

%% Example 2
Pdc = 0;
Vdc = 0;
%%
% When there is no DC power the inverter still draws power from the grid.
% This is expressed as a negative power output.
ACPower = pvl_adrinverter(Inverter, Vdc, Pdc)

%% References
%%
% * [1] A. Driesse, P. Jain, S. Harrison, Beyond the Curves: Modeling the Electrical Efficiency 
%       of Photovoltaic Inverters, 33rd IEEE PVSC, San Diego, CA 2008.
%% See Also 
% <pvl_snlinverter_help.html |pvl_snlinverter|>
%
%%
% Copyright 2014 Sandia National Laboratories

