function ACPower = pvl_adrinverter(Inverter, Vdc, Pdc)
% PVL_ADRINVERTER Converts DC power and voltage to AC power 
% using Anton Driesse's Grid-Connected PV Inverter efficiency model
% 
% Syntax
%   ACPower = pvl_adrinverter(Inverter, Vdc, Pdc)
%
% Description
%   Determine the AC power output of an inverter given the DC voltage, DC
%   power, and appropriate Model parameters. This function is similar to
%   the function pvl_snlinverter so that the two can be easily substituted
%   and compared, however the parameters required by the two are function
%   are very different.  (See the references for details.)
%   
%   The output, ACPower, is clipped at the maximum power output, and
%   gives a negative power during low-input power conditions, but does NOT
%   account for maximum power point tracking voltage windows nor maximum
%   current or voltage limits on the inverter.
%
% Inputs
%   Inverter - A struct defining the inverter to be used, giving the
%              inverter performance parameters according to the model
%              developed by Anton Driesse[1]. A set of inverter
%              performance parameters are provided with PV_LIB at
%              <PV_LIB Folder>/Required Data/<Current Driesse database name>
%
%   Required struct components are:
%   Inverter.Pacmax - The maximum AC output power value, used to clip the output if needed,(W) 
%   Inverter.Pnom   - The nominal power value used to normalize all power values, 
%                     typically the DC power needed to produce maximum AC power output,(W) 
%   Inverter.Vnom   - The nominal DC voltage value used to normalize DC voltage values,
%                     typically the level at which the highest efficiency is achieved, (V)
%   Inverter.Pnt    - ac-power consumed by inverter at night (night tare) to maintain 
%                     circuitry required to sense PV array voltage, (W)
%   Inverter.ADRCoefficients - This is a vector of 9 coefficients that capture the influence
%                   - of input voltage and power on inverter losses, and thereby efficiency
%
%   Vdc      - A scalar or vector of DC voltages, in volts, which are provided
%              as input to the inverter. If Vdc and Pdc are vectors, they must be 
%              of the same size. Vdc must be >= 0. (V)
%   Pdc      - A scalar or vector of DC powers, in watts, which are provided
%              as input to the inverter. If Vdc and Pdc are vectors, they must be 
%              of the same size. Pdc must be >= 0. (W)
%
% Outputs
%   ACPower - a column vector of modeled AC power output given the input 
%     DC voltage, Vdc, and input DC power, Pdc. When ACPower would be 
%     greater than Pmax, it is set to Pmax to represent inverter 
%     "clipping". When ACPower would be less than -Pnt (energy consumed rather
%     than produced) then ACPower is set to -Pnt to represent nightly 
%     power losses. ACPower is not adjusted for maximum power point
%     tracking (MPPT) voltage windows or maximum current limits of the
%     inverter.
%
% Reference:
%   [1] Beyond the Curves: Modeling the Electrical Efficiency 
%       of Photovoltaic Inverters, PVSC 2008, Anton Driesse et. al.
%
% See also
%   PVL_SAPM PVL_SNLINVERTER

p = inputParser;
p.addRequired('Inverter',@(x) isstruct(x))
p.addRequired('Vdc', @(x) all(isnumeric(x) & x>=0 & isvector(x)));
p.addRequired('Pdc', @(x) all(isnumeric(x) & x>=0 & isvector(x)));
p.parse(Inverter, Vdc, Pdc);

% Copy the model parameters from the inverter data structure for convenience
Pnom = p.Results.Inverter.Pnom;
Vnom = p.Results.Inverter.Vnom;
Pacmax = p.Results.Inverter.Pacmax;
Pnt  = p.Results.Inverter.Pnt;

b = p.Results.Inverter.ADRCoefficients;

% Make sure the inputs are columns
if isrow(Vdc) 
    Vdc = Vdc'; 
end;

if isrow(Pdc)
    Pdc = Pdc';
end

% Calculate the normalized or per unit dc voltage and power
vdc = Vdc / Vnom;
pdc = Pdc / Pnom;

% Define a function to expand the model polynomial
pv2X = @(p, v) [(p.^0) (p) (p.^2) (v-1) (p.*(v-1)) (p.^2.*(v-1)) (1./v-1) (p.*(1./v-1)) (p.^2.*(1./v-1))];

% Apply the model to calculate the normalized power loss
X = pv2X (pdc, vdc);
ploss = X * b';

% Calculate the output power scaled up to watts again
ACPower = Pnom * (pdc - ploss);

% Apply constraints 
Pnt  = -1 * abs(Pnt);           % Make Pnt a negative number (power loss)
ACPower((ACPower < Pnt) | (Vdc == 0))  = Pnt;  % Inverter night tare losses
ACPower(ACPower > Pacmax) = Pacmax; % Inverter clipping at maximum rated AC Power
