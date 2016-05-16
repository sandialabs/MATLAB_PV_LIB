function DCPower = pvl_huld(Module, Ee, Tm)
% PVL_HULD Calculates DC power using the Huld PV module model
% 
% Syntax
%   DCPower = pvl_huld(Module, G, Tm)
%
% Description
%   Determine the DC power output of PV module from in-plane effective
%   irradiance and cell temperature.
%
% Inputs
%   Module -  A struct defining the module being modeled and containing 
%             the following parameters:
%       Module.k - 1x6 vector with the coefficients k1 - k6
%       Module.Pmp0 - Power at STC reference condition (W)
%
%   Ee     - A scalar or vector of effective irradiance (suns)
%   Tm     - A scalar or vector of module temperature (C)
%
% Outputs
%   DCPower - a column vector of modeled DC power output
%
% Notes:
%   The Huld model [1] has the form of a polynomial motivated
%   by multiplying the separate equations for Vmp and Imp in the Sandia 
%   Array Performance Model [2].  The Huld model differs
%   from SAPM in that it uses Tm in place of Tc (cell temperature) 
%   The Huld model requires 7 parameters: Pmp0 (Pmp at STC) and 
%   six empirical coefficients k1 - k6 which can be estimated by:
%      - computation from coefficients for the SAPM
%      - fitting the model to efficiency vs. Ee and Tm data [1].
%   This implementation specifies Ee as input, as was done in [1] although 
%   not explicitly stated.
%
% Reference:
%   [1] A power-rating model for crystalline silicon PV modules, T. Huld,
%   G. Friesen, A. Skoczek, R. Kenny, T. Sample, M. Field, E. Dunlop, Solar
%   Energy Materials and Solar Cells 95(2011), pp 3359-3369
%
% See also
%   PVL_SAPM

p = inputParser;
p.addRequired('Module',@(x) isstruct(x))
p.addRequired('Ee', @(x) all(isnumeric(x) & x>=0 & isvector(x)));
p.addRequired('Tm', @(x) all(isnumeric(x) & x>=0 & isvector(x)));
p.parse(Module, Ee, Tm);

%   Define Constants
Const.T0 = 25;          %Reference temperature (25 deg C)

% Make sure the inputs are columns
Ee = Ee(:);
Tm = Tm(:);

% Calculate the DC power
DCPower = Ee.*(Module.Pmp0 + Module.k(1)*log(Ee) + Module.k(2)*log(Ee).^2 ...
    + Module.k(3)*(Tm - Const.T0) + Module.k(4)*(Tm - Const.T0).*log(Ee) ...
    + Module.k(5)*(Tm - Const.T0).*log(Ee).^2 + Module.k(6)*(Tm - Const.T0).^2);

