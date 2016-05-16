function [IL, I0, Rs, Rsh, nNsVth] = pvl_calcparams_PVsyst(Ee,Tcell,alpha_isc,ModuleParameters, varargin)
% PVL_CALCPARAMS_PVSYST Calculates five parameters for an IV curve using the PVsyst model.
%
% Syntax
%   [IL, I0, Rs, Rsh, nNsVth] = pvl_calcparams_PVsyst(S, Tcell, alpha_isc, ModuleParameters)
%   [IL, I0, Rs, Rsh, nNsVth] = pvl_calcparams_PVsyst(S, Tcell, alpha_isc, ModuleParameters)
%   [IL, I0, Rs, Rsh, nNsVth] = pvl_calcparams_PVsyst(S, Tcell, alpha_isc, ModuleParameters, Sref)
%   [IL, I0, Rs, Rsh, nNsVth] = pvl_calcparams_PVsyst(S, Tcell, alpha_isc, ModuleParameters, Sref, Tref)
%   [IL, I0, Rs, Rsh, nNsVth] = pvl_calcparams_PVsyst(S, Tcell, alpha_isc, ModuleParameters, 'Sref', Sref, 'Tref', Tref)
%   
% Description
%   Applies the temperature and irradiance corrections to calculate the
%   five parameters for an IV curve according to the PVsyst model [1, 2, 3].
%   The results of this procedure may be used to determine the IV curve
%   at effective irradiance = S, cell temperature = Tcell.
%
% Input Parameters:
%   Ee - The effective irradiance (in W/m^2) absorbed by the module. Ee must be >= 0.
%      May be a vector of irradiances, but must be the same size as all
%      other input vectors. Due to a division by Ee in the script, any Ee
%      value equal to 0 will be set to 1E-10.
%   Tcell - The average cell temperature of cells within a module in C.
%      Tcell must be >= -273.15. May be a vector of cell temperatures, but 
%      must be the same size as all other input vectors.
%   alpha_isc - The short-circuit current temperature coefficient of the 
%      module in units of A/C (or A/K).
%   ModuleParameters - a struct with parameters describing PV module
%     performance at reference conditions. The ModuleParameters struct 
%     must contain (at least) the following fields:
%      ModuleParameters.gamma_ref - diode (ideality) factor parameter at
%          reference conditions (unitless).
%      ModuleParameters.mugamma - temperature dependence of gamma (1/C)
%      ModuleParameters.IL_ref - Light-generated current (or photocurrent) 
%          in amperes at reference conditions. This value is referred to 
%          as Iphi in some literature.
%      ModuleParameters.I0_ref - diode reverse saturation current in amperes, 
%          under reference conditions.
%      ModuleParameters.Rsh_ref - shunt resistance under reference conditions (ohms)
%      ModuleParameters.Rsh0 - shunt resistance at zero irradiance (ohms)
%      ModuleParameters.Rshexp - exponential factor defining decrease in 
%          Rsh with increasing effective irradiance
%      ModuleParameters.Rs_ref - series resistance under reference conditions (ohms)
%      ModuleParameters.eG - The energy bandgap at reference temperature (in eV). 1.121 eV
%      for silicon. eG must be >0.
%   Sref - Optional reference effective irradiance in W/m^2. If omitted, a value of
%      1000 W/m^2 is used.
%   Tref - Optional reference cell temperature in C. If omitted, a value of
%      25 C is used.
%  
% Output:   
%   IL - Light-generated current in amperes at irradiance=S and 
%      cell temperature=Tcell. 
%   I0 - Diode saturation curent in amperes at irradiance S and cell temperature Tcell. 
%   Rs - Series resistance in ohms at irradiance S and cell temperature Tcell.
%   Rsh - Shunt resistance in ohms at irradiance S and cell temperature Tcell.
%   nNsVth - modified diode (ideality) factor at irradiance S and cell temperature
%      Tcell. nNsVth is the product of the usual diode (ideality) factor (gamma),
%      the number of series-connected cells in the module (Ns), and the thermal voltage 
%      of a cell in the module (Vth) at a cell temperature of Tcell.
%
%
% Sources:
%
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
% 
% See also
%   PVL_SINGLEDIODE    


p = inputParser;
p.addRequired('Ee',@(x) all(x>=0) & isnumeric(x) & isvector(x) );
p.addRequired('Tcell',@(x) all(x>=-273.15) & isnumeric(x) & isvector(x) );
p.addRequired('alpha_isc', @(x) (isnumeric(x) & isvector(x)));
p.addRequired('ModuleParameters', @(x) (isstruct(x)));
p.addOptional('Sref',1000, @(x) (all(x>0)& isnumeric(x) & isvector(x)));
p.addOptional('Tref',25, @(x) (all(x>-273.15) & isnumeric(x) & isvector(x)));
p.parse(Ee, Tcell, alpha_isc, ModuleParameters, varargin{:});

Ee = p.Results.Ee(:);
Sref = p.Results.Sref(:);
Tref = p.Results.Tref(:);
Tcell = p.Results.Tcell(:);
alpha_isc = p.Results.alpha_isc(:);
ModuleParameters = p.Results.ModuleParameters;

%Reference parameter should be taken and input from CEC database or similar
gamma_ref=ModuleParameters.gamma_ref(:);
mugamma=ModuleParameters.mugamma(:);
IL_ref=ModuleParameters.IL_ref(:);
I0_ref=ModuleParameters.I0_ref(:);
Rsh_ref=ModuleParameters.Rsh_ref(:);
Rsh0=ModuleParameters.Rsh0(:);
Rshexp=ModuleParameters.Rshexp(:);
Rs_ref=ModuleParameters.Rs_ref(:);
eG=ModuleParameters.eG(:);
% 
VectorSizes = [numel(Ee), numel(Tcell),...
    numel(Sref), numel(Tref), numel(IL_ref),...
    numel(I0_ref), numel(Rsh_ref), numel(Rsh0), numel(Rshexp),...
    numel(Rs_ref), numel(eG), numel(alpha_isc)];
MaxVectorSize = max(VectorSizes);
if not(all((VectorSizes==MaxVectorSize) | (VectorSizes==1)))
    error(['Input vectors S, Tcell, alpha_isc, Sref (if used), '...
        'Tref (if used), and all used components of ModuleParameters must '... 
        'either be scalars or vectors of the same length.']);
end

q=1.6021766E-19; % elementary charge in units of coulomb
k=1.3806488e-23; %Boltzman's constant in units of J/K
%k = 8.617332478e-5; % Boltzmann constant in units of eV/K


Tref_K=Tref+273.15; % Reference cell temperature in Kelvin
Tcell_K=Tcell+273.15; % cell temperature in Kelvin


gamma = gamma_ref+ mugamma*(Tcell - Tref); % Equation 8 in [1]
nNsVth = gamma*ModuleParameters.Ns*k/q.*Tcell_K;
IL=Ee./Sref.*(IL_ref+alpha_isc.*(Tcell_K-Tref_K)); 
IL(Ee <= 0) = 0; % If there is no light then no current is produced
% Equation 5 in [1]
I0=I0_ref.*((Tcell_K./Tref_K).^3).*exp((q*eG./(k.*gamma).*((1./Tref_K)-(1./Tcell_K))));
I0(IL==0) = 0; % If there is no light-generated current, there is no reverse saturation current
Rsh_base = max((Rsh_ref-Rsh0*exp(-Rshexp))./(1-exp(-Rshexp)),0); % Equation 7 in [1]
Rsh=Rsh_base+(Rsh0-Rsh_base).*exp(-Rshexp*(Ee./Sref)); % Equation 6 in [1]
Rsh(Ee <= 0) = inf; % Rsh is undefined if there is no current
Rs = Rs_ref;
