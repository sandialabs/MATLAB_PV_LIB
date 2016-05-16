function [IL, I0, Rs, Rsh, nNsVth] = pvl_calcparams_desoto(S,Tcell,alpha_isc,ModuleParameters, EgRef, dEgdT, varargin)
% PVL_CALCPARAMS_DESOTO Calculates five parameters to compute IV curves using the De Soto model.
%
% Syntax
%   [IL, I0, Rs, Rsh, nNsVth] = pvl_calcparams_desoto(S, Tcell, alpha_isc, ModuleParameters, EgRef, dEgdt)
%   [IL, I0, Rs, Rsh, nNsVth] = pvl_calcparams_desoto(S, Tcell, alpha_isc, ModuleParameters, EgRef, dEgdt)
%   [IL, I0, Rs, Rsh, nNsVth] = pvl_calcparams_desoto(S, Tcell, alpha_isc, ModuleParameters, EgRef, dEgdt, Sref)
%   [IL, I0, Rs, Rsh, nNsVth] = pvl_calcparams_desoto(S, Tcell, alpha_isc, ModuleParameters, EgRef, dEgdt, Sref, Tref)
%   [IL, I0, Rs, Rsh, nNsVth] = pvl_calcparams_desoto(S, Tcell, alpha_isc, ModuleParameters, EgRef, dEgdt, 'Sref', Sref, 'Tref', Tref)
%   
% Description
%   Applies the temperature and irradiance corrections to calculate the
%   five parameters for an IV curve according to the De Soto model [1].
%   The results of this procedure may be used to determine the IV curve
%   at effective irradiance S and cell temperature Tcell.
%
% Input Parameters:
%   S - The effective irradiance (in W/m^2) absorbed by the module. S must be >= 0.
%      May be a vector the same size as Tcell. Due to a division by S in the script, any S
%      value equal to 0 will be set to 1E-10.
%   Tcell - The average cell temperature of cells within a module in C.
%      Tcell must be >= -273.15. May be a vector of the same size as S.
%   alpha_isc - The short-circuit current temperature coefficient of the 
%      module in units of A/C (or A/K).
%   ModuleParameters - a struct with parameters describing PV module
%     performance at reference conditions according to DeSoto's paper. 
%     Parameters may be generated or found by lookup. For ease of use, a library 
%     of modules  A SAM library of modules is provided with PV_LIB (|\Required Data\CECModuleDatabaseSAM2014.1.14.mat|),
%     which may be read by the SAM library reader
%     function <pvl_SAMLibraryReader_CECModules_help.html |pvl_SAMLibraryReader_CECModules|>. may be found within the System Advisor Model (SAM) [2].
%     The SAM library has been provided as a .mat file, or may
%     be read and converted into a .mat file by the SAM library reader
%     functions. The ModuleParameters struct must contain (at least) the 
%     following fields:
%      ModuleParameters.a_ref - modified diode ideality factor parameter at
%          reference conditions, a_ref can be calculated from the
%          usual diode ideality factor (n), number of cells in series (Ns),
%          and cell temperature (Tcell) per equation (2) in [1].
%      ModuleParameters.IL_ref - Light-generated current (or photocurrent) 
%          in amperes at reference conditions. This value is referred to 
%          as Iphi in some literature.
%      ModuleParameters.I0_ref - diode reverse saturation current in amperes, 
%          under reference conditions.
%      ModuleParameters.Rsh_ref - shunt resistance under reference conditions (ohms)
%      ModuleParameters.Rs_ref - series resistance under reference conditions (ohms)
%   EgRef - The energy bandgap at reference temperature (in eV). 1.121 eV
%      for silicon. EgRef must be >0.
%   dEgdT - The temperature dependence of the energy bandgap at SRC (in 1/C).
%      May be either a scalar value (e.g. -0.0002677 as in [1]) or a
%      vector of dEgdT values corresponding to each input condition (this
%      may be useful if dEgdT is a function of temperature).
%   Sref - Optional reference irradiance in W/m^2. If omitted, a value of
%      1000 is used.
%   Tref - Optional reference cell temperature in C. If omitted, a value of
%      25 C is used.
%  
% Output:   
%   IL - Light-generated current in amperes at irradiance=S and 
%      cell temperature=Tcell. 
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
%    If the reference parameters in the ModuleParameters struct are read
%    from a database or library of parameters (e.g. System Advisor Model),
%    it is important to use the same EgRef and dEgdT values that
%    were used to generate the reference parameters, regardless of the 
%    actual bandgap characteristics of the semiconductor. For example, in 
%    the case of the System Advisor Model library, created as described in 
%    [3], EgRef and dEgdT for all modules were 1.121 and -0.0002677,
%    respectively.
%
% Sources:
%
% [1] W. De Soto et al., "Improvement and validation of a model for
%     photovoltaic array performance", Solar Energy, vol 80, pp. 78-88,
%     2006.
%
% [2] System Advisor Model web page. https://sam.nrel.gov.
%
% [3] A. Dobos, "An Improved Coefficient Calculator for the California
%     Energy Commission 6 Parameter Photovoltaic Module Model", Journal of
%     Solar Energy Engineering, vol 134, 2012.
%
% [4] O. Madelung, "Semiconductors: Data Handbook, 3rd ed." ISBN
%     3-540-40488-0
%
% See also
%   PVL_SAPM   PVL_SAPMCELLTEMP   PVL_SINGLEDIODE    
%      PVL_SAMLIBRARYREADER_CECMODULES


%  This table of reference bandgap energies (EgRef), bandgap energy
%  temperature dependence (dEgdT), and "typical" airmass response (M) is
%  provided purely as reference to those who may generate their own
%  reference module parameters (a_ref, IL_ref, I0_ref, etc.) based upon the
%  various PV semiconductors. Again, we stress the importance of
%  using identical EgRef and dEgdT when generation reference
%  parameters and modifying the reference parameters (for irradiance,
%  temperature, and airmass) per DeSoto's equations.
%  -----------------------------------------------------------------------
%
%   Silicon (Si):
%       EgRef = 1.121
%       dEgdT = -0.0002677
%       Source = Reference 1
%   Cadmium Telluride (CdTe):
%       EgRef = 1.475
%       dEgdT = -0.0003
%       Source = Reference 4
%   Copper Indium diSelenide (CIS):
%       EgRef = 1.010
%       dEgdT = -0.00011
%       Source = Reference 4
%   Copper Indium Gallium diSelenide (CIGS):
%       EgRef = 1.15
%       dEgdT = ????
%       Source = Wikipedia
%   Gallium Arsenide (GaAs):
%       EgRef = 1.424
%       dEgdT = -0.000433
%       Source = Reference 4

p = inputParser;
p.addRequired('S',@(x) all(x>=0) & isnumeric(x) & isvector(x) );
p.addRequired('Tcell',@(x) all(x>=-273.15) & isnumeric(x) & isvector(x) );
p.addRequired('alpha_isc', @(x) (isnumeric(x)));
p.addRequired('ModuleParameters', @(x) (isstruct(x)));
p.addRequired('EgRef', @(x) (isnumeric(x) & all(x>0)));
p.addRequired('dEgdT', @(x) (isnumeric(x)));
p.addOptional('Sref',1000, @(x) (all(x>0)& isnumeric(x)));
p.addOptional('Tref',25, @(x) (all(x>-273.15) & isnumeric(x)));
p.parse(S, Tcell, alpha_isc, ModuleParameters, EgRef, dEgdT, varargin{:});

S = p.Results.S(:);

Sref = p.Results.Sref(:);
Tref = p.Results.Tref(:);
Tcell = p.Results.Tcell(:);
alpha_isc = p.Results.alpha_isc(:);
ModuleParameters = p.Results.ModuleParameters;
dEgdT = p.Results.dEgdT(:);
EgRef = p.Results.EgRef(:);

a_ref=ModuleParameters.a_ref(:);
IL_ref=ModuleParameters.IL_ref(:);
I0_ref=ModuleParameters.I0_ref(:);
Rsh_ref=ModuleParameters.Rsh_ref(:);
Rs_ref=ModuleParameters.Rs_ref(:); 
 
VectorSizes = [numel(S), numel(Tcell)];
MaxVectorSize = max(VectorSizes);
if not(all((VectorSizes==MaxVectorSize) | (VectorSizes==1)))
    error(['Input vectors S and Tcell must '... 
        'either be scalars or vectors of the same length.']);
end

%k=1.3806488e-23; %Boltzman's constant in units of J/K
k = 8.617332478e-5; % Boltzmann constant in units of eV/K

Tref_K=Tref+273.15; % Reference cell temperature in Kelvin
Tcell_K=Tcell+273.15; % cell temperature in Kelvin

E_g=EgRef.*(1+dEgdT.*(Tcell_K-Tref_K)); % Equation 10 in [1]
nNsVth=a_ref.*(Tcell_K./Tref_K); % Equation 8 in [1]
IL=S./Sref.*(IL_ref+alpha_isc.*(Tcell_K-Tref_K)); 
IL(S <= 0) = 0; % If there is no light then no current is produced
I0=I0_ref.*((Tcell_K./Tref_K).^3).*exp((EgRef./(k.*Tref_K))-(E_g./(k.*Tcell_K)));
I0(IL==0) = 0; % If there is no light-generated current, there is no reverse saturation current
Rsh=Rsh_ref.*(Sref./S); % Equation 12 in [1]
Rsh(S <= 0) = inf; % Rsh is undefined if there is no current
Rs = Rs_ref;
