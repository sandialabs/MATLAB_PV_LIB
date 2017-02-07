function [M] = pvl_FSspeccorr(Pwat,AMa,varargin)
% pvl_FSspeccorr  Spectral mismatch modifier based on precipitable water and
%  absolute (pressure corrected) airmass.
%
% Syntax:
%   [M] = pvl_FSspeccorr(Pwat, AMa, pvModType)
%   [M] = pvl_FSspeccorr(Pwat, AMa, custCoeff)
%
% Description:
%
%   Estimates a spectral mismatch modifier M representing the effect on
%   module short circuit current of variation in the spectral irradiance.
%   M is estimated from absolute (pressure currected) air mass, AMa, and
%   precipitable water, Pwat, using the following function:
%
%   M = coeff(1) + coeff(2)*AMa  + coeff(3)*Pwat  + coeff(4)*AMa.^.5
%           + coeff(5)*Pwat.^.5 + coeff(6)*AMa./Pwat.^0.5         (1)
%
%   Default coefficients are determined for several cell types with
%   known quantum efficiency curves, by using the Simple Model of the
%   Atmospheric Radiative Transfer of Sunshine (SMARTS) [1].
%   Using SMARTS, spectrums are simulated with all combinations of AMa
%   and Pwat where:
%       *   0.1 cm <= Pwat <= 5 cm
%       *   1.0 <= AMa <= 5
%       *   Spectral range is limited to that of CMP11 (280 nm to 2800 nm)
%       *   spectrum simulated on a plane normal to the sun
%       *   All other parameters fixed at G173 standard
%   From these simulated spectra, M is calculated using the known quantum
%   efficiency curves. Multiple linear regression is then applied to fit
%   Eq. 1 to determine the coefficients for each module.
%
%  Function pvl_FSspeccorr was developed by Mitchell Lee and Alex Panchula,
%   at First Solar, 2015. Detailed description of the spectral correction
%   can be found in [2]. Additional validation and testing of the model
%   can be found in [3].
%
% Inputs:
%   Pwat - atmospheric precipitable water (cm). Can be
%       entered as a vector.
%   AMa - absolute (pressure corrected) airmass, as a vector of the same
%       length as Pwat
%   pvModType - a string specifying a cell type. Can be lower or upper case
%       letters.  Admits values of 'cdte', 'monosi'='xsi', 'multisi'='polysi'.
%       If provided, this input
%       selects coefficients for the following default modules:
%           'cdte' - coefficients for First Solar Series 4-2 CdTe modules.
%           'monosi','xsi' - coefficients for First Solar TetraSun modules.
%           'multisi','polysi' - coefficients for multi-crystalline silicon
%               modules. The module used to calculate the spectral
%               correction coefficients corresponds to the Mult-crystalline
%               silicon Manufacturer 2 Model C from [4].
%           'cigs' - coefficients for anonymous copper indium gallium selenide
%               PV module. Lower and upper limits of QE are 350 nm and 1300 nm,
%               respectively. Please note that the QE of CIGS modules
%               can vary significantly depending on the PV manufacture and
%               vintage. Spectral Response of module module used to derive
%               CIGS coefficients can be found in [3]. 
%           'asi' - coefficients for anonymous amorphous silicon PV module.
%               Lower and upper limits of QE are 280 nm and 800 nm,
%               respectively. Please note that the QE of a-Si modules
%               can vary significantly depending on the PV manufacture and
%               vintage.  Spectral Response of module module used to derive
%               a-Si coefficients can be found in [3]. 
%   custCoeff - allows for entry of user defined spectral correction
%       coefficients. Coefficients must be entered as a numeric row or
%       column vector of length 6. Derivation of coefficients requires use
%       of SMARTS and PV module quantum efficiency curve. Useful for modeling
%       PV module types which are not included as defaults, or to fine tune
%       the spectral correction to a particular mono-Si, multi-Si, or CdTe
%       PV module. Note that the parameters for modules with very
%       similar QE should be similar, in most cases limiting the need for
%       module specific coefficients.
%
%
% Output:
%   M - spectral mismatch factor (unitless) which is can be multiplied
%       with broadband irradiance reaching a module's cells to estimate
%       effective irradiance, i.e., the irradiance that is converted
%       to electrical current.
%
% References:
% [1]   Gueymard, Christian. SMARTS2: a simple model of the atmospheric
%           radiative transfer of sunshine: algorithms and performance
%           assessment. Cocoa, FL: Florida Solar Energy Center, 1995.
% [2]   Lee, Mitchell, and Panchula, Alex. "Spectral Correction for
%           Photovoltaic Module Performance Based on Air Mass and Precipitable
%           Water." IEEE Photovoltaic Specialists Conference, Portland, 2016
% [3]   Schweiger, M. and Hermann, W, Influence of Spectral Effects on 
%           Energy Yield of Different PV Modules: Comparison of Pwat and 
%           MMF Approach, TUV Rheinland Energy GmbH report 21237296.003, January 2017
% [4]   Marion, William F., et al. User's Manual for Data for Validating
%           Models for PV Module Performance. National Renewable Energy Laboratory, 2014.
%           http://www.nrel.gov/docs/fy14osti/61610.pdf


% Correct for AMa and Pwat having transposed dimensions
if isrow(AMa)
    AMa = AMa';
end

if isrow(Pwat)
    Pwat = Pwat';
end

% --- Screen Input Data ---

% *** Pwat ***
% Replace Pwat Values below 0.1 cm with 0.1 cm to prevent model from
% diverging
if min(Pwat) < 0.1
    Pwat(Pwat < 0.1) = 0.1;
    warning(['Exceptionally low Pwat values replaced with 0.1 cm to prevent',...
        ' model divergence']);
end

% Warn user about Pwat data that is exceptionally high
if max(Pwat) > 8
    warning(['Exceptionally high Pwat values. Check input data:', ...
        ' model may diverge in this range']);
end

% *** AMa ***
% Replace Extremely High AM with AM 10 to prevent model divergence
% AM > 10 will only occur very close to sunset
if max(AMa) > 10
    AMa(AMa > 10) = 10;
end

% Warn user about AMa data that is exceptionally low
if min(AMa) < 0.58
    warning(['Exceptionally low air mass: ',...
        'model not intended for extra-terrestrial use'])
    % pvl_absoluteairmass(1,pvl_alt2pres(4340)) = 0.58
    % Elevation of Mina Pirquita, Argentian = 4340 m. Highest elevation city
    % with population over 50,000.
end

% If user input is a character array, use appropriate default coefficients.
if ischar(varargin{1})
    modType = lower(varargin{1});
    modType = regexprep(modType,'[^a-zA-Z]','');
    switch modType
        case 'cdte'
            % Coefficients for First Solar Series 4-2 (and later) modules.
            % For modeling the performance of earlier CdTe module series,
            % use the coefficients that are commented out
            % [0.79418,-0.049883,-0.013402,0.16766,0.083377,-0.0044007];
            coeff = [0.86273, -0.038948, -0.012506, 0.098871, 0.084658, -0.0042948];
        case {'monosi','xsi'}
            % Coefficients for First Solar TetraSun Modules
            coeff = [0.85914, -0.020880, -0.0058853, 0.12029, 0.026814, -0.0017810];
        case {'polysi','multisi'}
            % Coefficients for Multi-Si: Manufacturer 2 Model C
            coeff = [0.84090, -0.027539, -0.0079224, 0.13570, 0.038024, -0.0021218];
        case 'cigs'
            coeff = [0.85252, -0.022314, -0.0047216, 0.13666, 0.013342, -0.0008945];
        case 'asi'
            coeff = [1.12094, -0.047620, -0.0083627, -0.10443, 0.098382,-0.0033818];
        otherwise
            error('Incorrect module type for use of default parameters')
    end
% User input coefficients
else
    coeff = varargin{1};
end

% Evaluate Spectral Shift
M = coeff(1) + coeff(2)*AMa  + coeff(3)*Pwat  + coeff(4)*AMa.^.5  + coeff(5)*Pwat.^.5 + coeff(6)*AMa./Pwat.^0.5;
end

