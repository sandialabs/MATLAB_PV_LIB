function Result = pvl_sapm(Module, Ee, celltemp)
% PVL_SAPM Performs Sandia PV Array Performance Model to get 5 points on IV curve given SAPM module parameters, Ee, and cell temperature
%
% Syntax
%   Result = pvl_sapm(Module, Ee, celltemp)
%
% Description
%   The Sandia PV Array Performance Model (SAPM) generates 5 points on a PV
%   module's I-V curve (Voc, Isc, Ix, Ixx, Vmp/Imp) according to
%   SAND2004-3535. Assumes a reference cell temperature of 25 C.
% Inputs
%   Module - a structure defining the SAPM performance parameters (see
%   below)
%   Ee - The effective irradiance incident upon the module (suns). Any Ee<0
%   are set to 0.
%   celltemp - The cell temperature (degrees C)
% Outputs
%   Result - A structure with Result.Isc, Result.Imp, Result.Ix,
%   Result.Ixx, Result.Voc, Result.Vmp, and Result.Pmp. Structure
%   components are vectors of the same size as Ee.
%
% Note - The particular coefficients from SAPM which are required in Module
% are:
%   Module.c - 1x8 vector with the C coefficients Module.c(1) = C0
%   Module.Isc0 - Short circuit current at reference condition (amps)
%   Module.Imp0 - Maximum power current at reference condition (amps)
%   Module.AlphaIsc - Short circuit current temperature coefficient at
%     reference condition (1/C)
%   Module.AlphaImp - Maximum power current temperature coefficient at
%     reference condition (1/C)
%   Module.BetaVoc - Open circuit voltage temperature coefficient at
%     reference condition (V/C)
%   Module.mBetaVoc - Coefficient providing the irradiance dependence for
%     the BetaVoc temperature coefficient at reference irradiance (V/C)
%   Module.BetaVmp - Maximum power voltage temperature coefficient at
%     reference condition
%   Module.mBetaVmp - Coefficient providing the irradiance dependence for
%     the BetaVmp temperature coefficient at reference irradiance (V/C)
%   Module.n - Empirically determined "diode factor" (dimensionless)
%   Module.Ns - Number of cells in series in a module's cell string(s)
%
% References
%   [1] King, D. et al, 2004, "Sandia Photovoltaic Array Performance Model", SAND Report
%   3535, Sandia National Laboratories, Albuquerque, NM
%
% See also PVL_SAPMMODULEDB PVL_SAPMCELLTEMP 

p = inputParser;
p.addRequired('Module',@(x) all(isstruct(x)));
p.addRequired('Ee', @(x) all(isvector(x) & (isnumeric(x) | isnan(x))));
p.addRequired('celltemp', @(x) all(isvector(x) & isnumeric(x)));
p.parse(Module,Ee,celltemp);


% Make sure inputs are of the same size
if size(Ee) ~= size(celltemp)
    error('Error in pvl_sapm. celltemp and Ee must be vectors of the same size')
end

%   Define Constants
Const.T0 = 25;          %Reference temperature (25 deg C)
Const.q = 1.60218E-19;  %Elementary charge (1.60218E-19 coulombs
Const.k = 1.38066E-23;  %Boltzmann's constant (1.38066E-23 J/K)

% Create output variables, set them to 0
Result.Isc = 0*Ee;
Result.Imp = 0*Ee;
Result.Voc = 0*Ee;
Result.Vmp = 0*Ee;
Result.Ix = 0*Ee;
Result.Ixx = 0*Ee;
Result.Pmp = 0*Ee;

% Create other variables, set them to 0
delta = 0*Ee;
BetaVoc = 0*Ee;
BetaVmp = 0*Ee;

Ee(Ee<0) = 0;
filter = (Ee >= 1E-3); % Don't perform SAPM on Ee values < 1E-3

Result.Isc(filter) = Module.Isc0 .* Ee(filter) .* (1+Module.AlphaIsc .*(celltemp(filter) - Const.T0));
Result.Imp(filter) = Module.Imp0.*(Module.c(1).* Ee(filter) + Module.c(2)*(Ee(filter).^2)).*(1 + Module.AlphaImp .* (celltemp(filter) - Const.T0));
BetaVoc(filter) = Module.BetaVoc + Module.mBetaVoc .* (1-Ee(filter));
delta(filter) = Module.n .* Const.k .* (celltemp(filter) + 273.15) ./ Const.q;
Result.Voc(filter) = (Module.Voc0 + Module.Ns .* delta(filter) .* log(Ee(filter)) + BetaVoc(filter) .* (celltemp(filter) - Const.T0));
BetaVmp(filter) = Module.BetaVmp + Module.mBetaVmp .* (1-Ee(filter));
Result.Vmp(filter) = (Module.Vmp0 + Module.c(3) .* Module.Ns .* delta(filter) .* log(Ee(filter))+ Module.c(4) .* Module.Ns .* (delta(filter) .* log(Ee(filter))).^2 + BetaVmp(filter) .* (celltemp(filter) - Const.T0));
Result.Pmp = Result.Imp .* Result.Vmp;
Result.Ix(filter) = Module.Ix0*(Module.c(5) .* Ee(filter) + Module.c(6) .* (Ee(filter)).^2) .*...
    (1 + Module.AlphaIsc.*(celltemp(filter)-Const.T0));
Result.Ixx(filter) = Module.Ixx0*(Module.c(7) .* Ee(filter) + Module.c(8) .* (Ee(filter)).^2) .*...
    (1 + Module.AlphaIsc.*(celltemp(filter)-Const.T0));
