function ACPower = pvl_snlinverter(Inverter, Vdc, Pdc)
% PVL_SNLINVERTER Converts DC power and voltage to AC power using Sandia's Grid-Connected PV Inverter model
% 
% Syntax
%   ACPower = pvl_snlinverter(Inverter, Vdc, Pdc)
%
% Description
%   Determine the AC power output of an inverter given the DC voltage, DC
%   power, and appropriate Sandia Grid-Connected Photovoltaic Inverter
%   Model parameters. The output, ACPower, is clipped at the maximum power
%   output, and gives a negative power during low-input power conditions,
%   but does NOT account for maximum power point tracking voltage windows
%   nor maximum current or voltage limits on the inverter. 
%
% Inputs
%   Inverter - A struct defining the inverter to be used, giving the
%     inverter performance parameters according to the Sandia
%     Grid-Connected Photovoltaic Inverter Model (SAND 2007-5036) [1]. A set of
%     inverter performance parameters are provided with PV_LIB at
%     <PV_LIB Folder>/Required Data/<Current SNLInverter database name>, 
%     or may be generated from a System Advisor Model (SAM) [2] library 
%     using the SAM library reader functions. 
%   
%   Required struct components are:
%   Inverter.Pac0 - AC-power output from inverter based on input power and voltage, (W) 
%   Inverter.Pdc0 - DC-power input to inverter, typically assumed to be equal to the PV array maximum
%               power, (W)
%   Inverter.Vdc0 - DC-voltage level at which the AC-power rating is achieved at the reference operating
%               condition, (V)
%   Inverter.Ps0 - DC-power required to start the inversion process, or self-consumption by inverter,
%               strongly influences inverter efficiency at low power levels, (W)
%   Inverter.C0 - Parameter defining the curvature (parabolic) of the relationship between ac-power and
%               dc-power at the reference operating condition, default value of zero gives a linear
%               relationship, (1/W)
%   Inverter.C1 - Empirical coefficient allowing Pdco to vary linearly with dc-voltage input, default value
%               is zero, (1/V)
%   Inverter.C2 - empirical coefficient allowing Pso to vary linearly with dc-voltage input, default value
%               is zero, (1/V)
%   Inverter.C3 - empirical coefficient allowing Co to vary linearly with dc-voltage input, default value is
%               zero, (1/V)
%   Inverter.Pnt - ac-power consumed by inverter at night (night tare) to maintain circuitry required to
%               sense PV array voltage, (W)
%   Vdc - A scalar or vector of DC voltages, in volts, which are provided
%     as input to the inverter. If Vdc and Pdc are vectors, they must be 
%     of the same size. Vdc must be >= 0.
%   Pdc - A scalar or vector of DC powers, in watts, which are provided
%     as input to the inverter. If Vdc and Pdc are vectors, they must be 
%     of the same size. Pdc must be >= 0.
%
% Outputs
%   ACPower - a column vector of modeled AC power output given the input 
%     DC voltage, Vdc, and input DC power, Pdc. When ACPower would be 
%     greater than Pac0, it is set to Pac0 to represent inverter 
%     "clipping". When ACPower would be less than Ps0 (startup power
%     required), then ACPower is set to -1*abs(Pnt) to represent nightly 
%     power losses. ACPower is not adjusted for maximum power point
%     tracking (MPPT) voltage windows or maximum current limits of the
%     inverter.
%
% Reference:
%   [1] (SAND2007-5036, "Performance Model for Grid-Connected Photovoltaic 
%   Inverters by D. King, S. Gonzalez, G. Galbraith, W. Boyson)
%
%   [2] System Advisor Model web page. https://sam.nrel.gov.
%
% See also
%   PVL_SAPM    PVL_SAMLIBRARYREADER_SNLINVERTERS   PVL_ADRINVERTER

p = inputParser;
p.addRequired('Inverter',@(x) isstruct(x))
p.addRequired('Vdc', @(x) all(isnumeric(x) & x>=0 & isvector(x)));
p.addRequired('Pdc', @(x) all(isnumeric(x) & x>=0 & isvector(x)));
p.parse(Inverter, Vdc, Pdc);

Pac0 = p.Results.Inverter.Pac0;
Pdc0 = p.Results.Inverter.Pdc0;
Vdc0 = p.Results.Inverter.Vdc0;
Ps0 = p.Results.Inverter.Ps0;
C0 = p.Results.Inverter.C0;
C1 = p.Results.Inverter.C1;
C2 = p.Results.Inverter.C2;
C3 = p.Results.Inverter.C3;
Pnt = p.Results.Inverter.Pnt;


A = Pdc0 .* (1 + C1 .* (Vdc - Vdc0));
B = Ps0 .* (1+ C2 .* (Vdc - Vdc0));
C = C0 .* (1 + C3 .* (Vdc - Vdc0));

ACPower = ((Pac0 ./ (A - B)) - C .* (A-B)) .* (Pdc - B)...
    + C .* (Pdc-B).^2;
    
ACPower(ACPower>Pac0) = Pac0; % Inverter clipping at maximum rated AC Power
ACPower(ACPower<Ps0)= -1.*abs(Pnt); % Inverter night tare losses
ACPower = ACPower(:);