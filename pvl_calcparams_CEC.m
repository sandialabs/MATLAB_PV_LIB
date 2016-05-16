function [IL, I0, Rs, Rsh, nNsVth] = pvl_calcparams_CEC(S,Tcell,ModuleParameters)
% PVL_CALCPARAMS_CEC Calculates five parameters for an IV curve using the CEC model.
%
% Syntax
%   [IL, I0, Rs, Rsh, nNsVth] = pvl_calcparams_CEC(S, Tcell, ModuleParameters)
%   
% Description
%   Applies the temperature and irradiance corrections to calculate the
%   five parameters for an IV curve according to the CEC model [1].
%   The results of this procedure may be used to determine the IV curve
%   at effective irradiance = S, cell temperature = Tcell.
%
% Input Parameters:
%   S - The effective irradiance (in W/m^2) absorbed by the module. S must be >= 0.
%      May be a vector of the same size as Tcell. Due to a division by S in the script, any S
%      value equal to 0 will be set to 1E-10.
%   Tcell - The average cell temperature of cells within a module in C.
%      Tcell must be >= -273.15. May be a vector of the same size as S.
%   ModuleParameters - a struct with parameters for the module. A library 
%     of parameters may be found within the System Advisor Model (SAM) [2].
%     The SAM library has been provided as a .mat file, or may
%     be read and converted into a .mat file by the SAM library reader
%     functions. The ModuleParameters struct must contain (at least) the 
%     following fields:
%       ModuleParameters.a_ref - modified diode ideality factor parameter at
%          reference conditions (units of V), a_ref can be calculated as 
%          a_ref = n Ns Vth, where n is the usual diode ideality factor (n),
%             Ns is the number of cells in series, and Vth is the thermal
%             voltage at STC cell temperature 298.15K.
%       ModuleParameters.IL_ref - Light-generated current (or photocurrent) 
%          in amperes at reference conditions. 
%       ModuleParameters.I0_ref - diode reverse saturation current in amperes, 
%          under reference conditions.
%       ModuleParameters.Rsh_ref - shunt resistance under reference conditions (ohms)
%       ModuleParameters.Rs_ref - series resistance under reference conditions (ohms)
%       ModuleParameters.adjust - an adjustment factor applied to the
%          reference value for the temperature coefficient for short circuit
%          current (percent)
%       ModuleParameters.alpha_sc - temperature coefficient for
%          short-circuit current at reference conditions (A/C)
%       ModuleParameters.Ns - number of cells in series (unitless)
%  
% Output:   
%   IL - Light-generated current in amperes at irradiance S and 
%      cell temperature Tcell. 
%   I0 - Diode saturation curent in amperes at irradiance S and cell temperature Tcell. 
%   Rs - Series resistance in ohms at irradiance S and cell temperature Tcell.
%   Rsh - Shunt resistance in ohms at irradiance S and cell temperature Tcell.
%   nNsVth - modified diode ideality factor at irradiance S and cell temperature
%      Tcell. Note that in source [1] nNsVth = a (equation 2). nNsVth is the 
%      product of the usual diode ideality factor (n), the number of 
%      series-connected cells in the module (Ns), and the thermal voltage 
%      of a cell in the module (Vth) at a cell temperature of Tcell.
%
% Notes:
%    In the case of the CEC model and the parameters in the System Advisor 
%    Model library, created as described in [3], EgRef and dEgdT for all 
%    modules are 1.121 and -0.0002677, respectively.
%
% Sources:
%
% [1] P. Gilman, SAM Photovoltaic Model Technical Reference, National
% Renewable Energy Laboratory (NREL) Technical Report NREL/TP-6A20-64102,
% May 2015
%
% [2] System Advisor Model web page. https://sam.nrel.gov.
%
% [3] A. Dobos, "An Improved Coefficient Calculator for the California
%     Energy Commission 6 Parameter Photovoltaic Module Model", Journal of
%     Solar Energy Engineering, vol 134, 2012.
%
%
% See also
%   PVL_CALCPARAMS_DESOTO  PVL_SINGLEDIODE    
%      PVL_SAMLIBRARYREADER_CECMODULES



p = inputParser;
p.addRequired('S',@(x) all(x>=0) & isnumeric(x) & isvector(x) );
p.addRequired('Tcell',@(x) all(x>=-273.15) & isnumeric(x) & isvector(x) );
p.addRequired('ModuleParameters', @(x) (isstruct(x)));
p.parse(S, Tcell, ModuleParameters);

S = p.Results.S(:);
Tcell = p.Results.Tcell(:);
ModuleParameters = p.Results.ModuleParameters;

a_ref=ModuleParameters.a_ref(:);
IL_ref=ModuleParameters.IL_ref(:);
I0_ref=ModuleParameters.I0_ref(:);
Rsh_ref=ModuleParameters.Rsh_ref(:);
Rs_ref=ModuleParameters.Rs_ref(:); 
adjust=ModuleParameters.adjust(:);
alpha_sc=ModuleParameters.alpha_sc(:);

VectorSizes = [numel(S), numel(Tcell), numel(a_ref), numel(IL_ref),...
    numel(I0_ref), numel(Rsh_ref), numel(Rs_ref), numel(alpha_sc)];
MaxVectorSize = max(VectorSizes);
if not(all((VectorSizes==MaxVectorSize) | (VectorSizes==1)))
    error(['Input vectors S, Tcell, and all used components of ModuleParameters must '... 
        'either be scalars or vectors of the same length.']);
end

%k=1.3806488e-23; %Boltzman's constant in units of J/K
k = 8.617332478e-5; % Boltzmann constant in units of eV/K

Sref = 1000;  % Reference effective irradiance in W/m2
Tref_K=25+273.15; % Reference cell temperature in Kelvin
Tcell_K=Tcell+273.15; % cell temperature in Kelvin
EgRef = 1.121;   % in eV, see [1]
dEgdT = -0.00002677;  % see [1]

% adjust temperature coefficient for short-circuit current, Eq. 9.11 of [1]
aIsc = alpha_sc .* ( 1 - adjust/100);

%These parameters (a, I_L, I_o, M, and R_sh) need to be vectors of equal
%length to the number of conditions (number of S,Tcell pairs).
E_g=EgRef.*(1+dEgdT.*(Tcell_K-Tref_K)); % Equation 10 in [1]
nNsVth=a_ref.*(Tcell_K./Tref_K); % Equation 8 in [1]
IL=S./Sref.*(IL_ref+aIsc.*(Tcell_K-Tref_K)); 
IL(S <= 0) = 0; % If there is no light then no current is produced
I0=I0_ref.*((Tcell_K./Tref_K).^3).*exp((EgRef./(k.*Tref_K))-(E_g./(k.*Tcell_K)));
I0(IL==0) = 0; % If there is no light-generated current, there is no reverse saturation current
Rsh=Rsh_ref.*(Sref./S); % Equation 12 in [1]
Rsh(S <= 0) = inf; % Rsh is undefined if there is no current
Rs = Rs_ref;
