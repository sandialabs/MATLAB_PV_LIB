%% pvl_snlinverter
% Sandia photovoltaic inverter model.
%
%% Syntax
%
% |ACPower = pvl_snlinverter(Inverter, Vdc, Pdc)| 
%% Description
% Computes the AC power output using the Sandia photovoltaic inverter model [1].
%
%% Inputs
%%
% * *|Vdc|* - a vector of modeled DC voltage.
% * *|Pdc|* - a vector of modeled DC power.
% * *|Inverter|* - A struct of parameters defining the inverter to be used. A library of
%     inverter performance parameters are provided with PV_LIB (|\Required 
%     Data\SandiaInverterDatabaseSAM2012.11.30.mat|), or a library may be
%     generated from a System Advisor Model (SAM) [2] library using the SAM
%     library reader function <pvl_SAMLibraryReader_SNLInverters_help.html |pvl_SAMLibraryReader_SNLInverters|>. 
%     Required struct components are:
%%
% * *|Inverter.Pac0|* - AC-power output from inverter based on input power
% and voltage, (W).
% * *|Inverter.Pdc0|* - DC-power input to inverter, typically assumed to be equal to the PV array maximum
% power, (W).
% * *|Inverter.Vdc0|* - DC-voltage level at which the AC-power rating is achieved at the reference operating
% condition, (V).
% * *|Inverter.Ps0|* - DC-power required to start the inversion process, or self-consumption by inverter,
% strongly influences inverter efficiency at low power levels, (W).
% * *|Inverter.C0|* - Parameter defining the curvature (parabolic) of the relationship between ac-power and
% dc-power at the reference operating condition, default value of zero gives a linear
% relationship, (1/W).
% * *|Inverter.C1|* - Empirical coefficient allowing Pdco to vary linearly with dc-voltage input, default value
% is zero, (1/V).
% * *|Inverter.C2|* - empirical coefficient allowing Pso to vary linearly with dc-voltage input, default value
% is zero, (1/V).
% * *|Inverter.C3|* - empirical coefficient allowing Co to vary linearly with dc-voltage input, default value is
% zero, (1/V).
% * *|Inverter.Pnt|* - ac-power consumed by inverter at night (night tare) to maintain circuitry required to
% sense PV array voltage, (W).
%
%% Outputs
%%
% * *|ACPower|* - a vector of modeled AC power output 
%%
% Note: When |ACPower| would be greater than |Pac0|, it is set to |Pac0| to represent
% inverter "clipping". When |ACPower| would be less than |Ps0| (startup power
% required), then |ACPower| is set to -1 * abs(|Pnt|) to represent tare 
% losses.

%% Example 1
load 'SandiaInverterDatabaseSAM2014.1.14.mat';
% PV Powered PVP2500 is entry #793
Inverter = SNLInverterDB(793)
%%
% Assume DC power is 1,000 W and voltage is 400 V
%
Pdc = 1000; %DC power is 1000 W
Vdc = 400; % DC violtage is 450 V
ACPower = pvl_snlinverter(Inverter, Vdc, Pdc)
%%
% Inverter efficiency can be calculated as:
Inverter_Efficiency = ACPower/Pdc

%% Example 2
Pdc = 0;
Vdc = 0;
%%
% When there is no DC power the inverter still draws power from the grid.
% This is expressed as a negative power output.
ACPower = pvl_snlinverter(Inverter, Vdc, Pdc)

%% References
%
% [1] King, D. et al., 2007. Performance Model for Grid-Connected Photovoltaic Inverters,
% SAND2007-5036, Sandia National Laboratories, Albuquerque, NM.
% <http://energy.sandia.gov/wp/wp-content/gallery/uploads/Performance-Model-for-Grid-Connected-Photovoltaic-Inverters.pdf
% Web Link>
%
% [2] System Advisor Model web page. <https://sam.nrel.gov Web Link>.

%% See Also 
% <pvl_sapm_help.html |pvl_sapm|>  
%%
% Copyright 2014 Sandia National Laboratories

