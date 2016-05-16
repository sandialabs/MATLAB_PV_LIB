
function [DNI, DHI, Kt] = pvl_reindl_1(GHI, Z, doy)
% PVL_REINDL_1 Estimate DNI and DHI from GHI using the Reindl_1 model
%
% Syntax
% [DNI, DHI, Kt] = pvl_reindl_1(GHI, Z, doy)
%
%
% Description
%   The Reindl-1 model estimates the diffuse fraction DF from global horizontal
%   irradiance through an empirical relationship between DF and the ratio of GHI to 
%   extraterrestrial irradiance, Kt.  pvl_reindl_1 uses the diffuse 
%   fraction to compute DHI. DNI is then estimated as DNI = (GHI - DHI)/cos(Z).
%
% Inputs:   
%   GHI - a scalar or vector of global horizontal irradiance in W/m^2. If GHI
%     is a vector it must be of the same size as all other vector inputs. 
%     GHI must be >=0.
%   Z - a scalar or vector of true (not refraction-corrected) zenith
%     angles in decimal degrees. If Z is a vector it must be of the
%     same size as all other vector inputs. Z must be >=0 and <=180.
%   doy - a scalar or vector of values providing the day of the year. If
%     doy is a vector it must be of the same size as all other vector inputs.
%     doy must be >= 1 and < 367.
%
% Output:   
%   DNI - the modeled direct normal irradiance in W/m^2.
%   DHI - the modeled diffuse horizontal irradiance in W/m^2.
%   Kt - Ratio of global to extraterrestrial irradiance on a horizontal
%     plane.
%
% Sources
%
% [1] Reindl DT, Beckman WA, Duffie JA. Diffuse fraction correlations.
% Solar Energy 1990;45:1-7
%   
%  
%
%
% See also PVL_DATE2DOY PVL_EPHEMERIS PVL_ALT2PRES PVL_DIRINT PVL_LOUCHE
% PVL_ORGILL_HOLLANDS PVL_REINDL_2 PVL_ERBS PVL_DISC

p = inputParser;
p.addRequired('GHI', @(x) all(isnumeric(x) & isvector(x) ));
p.addRequired('Z', @(x) (all(isnumeric(x) & x<=180 & x>=0 & isvector(x))));
p.addRequired('doy', @(x) (all(isnumeric(x) & isvector(x) & x>=1 & x<367)));
p.parse(GHI,Z,doy);


% Initialize variables
GHI = GHI.*ones(max([numel(GHI) numel(Z) numel(doy)]),1);
Z=Z(:);
doy=doy(:);
DF = zeros(length(GHI),1);

% The following code and comments utilize the model's calculations for
% extraterrestrial radiation. I'm not exactly sure what Maxwell was using
% as the "Eccentricity of the Earth's orbit", but I can't figure out what
% to put in to make the equations work out correctly.
% % It is unclear in Maxwell whether the trigonometric functions should
% % operate on degrees or radians. Spencer's work also does not explicitly
% % state the units to determine re (denoted as 1/r^2 in Spencer's work).
% % However, Spencer uses radian measures for earlier calculations, and it is
% % assumed to be similar for this calculation. In either case (radians or
% % degrees) the difference between the two methods is approximately 0.0015%.
% re = 1.00011 + 0.034221 .* cos(Eccentricity) + (0.00128) .* sin(Eccentricity)...
%     +0.000719.*cos(2.*Eccentricity) + (7.7E-5).*sin(2.*Eccentricity);
% I0= re.*Hextra;

HExtra = pvl_extraradiation(doy);
I0h= HExtra.*cosd(Z);

Kt = GHI./(I0h); % This Z needs to be the true Zenith angle, not apparent (to get extraterrestrial horizontal radiation)
Kt(Kt<0) = 0;


for i = 1:length(GHI)
        % For Kt <= 0.30, set the diffuse fraction
    if Kt(i)<=.30
        DF(i) = 1.02 - 0.248*Kt(i);
        
        % For Kt > 0.30 and Kt <= 0.78, set the diffuse fraction
    elseif Kt(i)>.30 && Kt(i)<.78
        DF(i) = 1.45 - 1.67*Kt(i);
        
        % For Kt > 0.78, set the diffuse fraction
    else
        DF(i) = 0.147;
    end
end


DHI = DF.*GHI;

DNI = (GHI - DHI)./(cosd(Z));

























