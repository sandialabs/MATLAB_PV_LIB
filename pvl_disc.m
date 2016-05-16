function [DNI, Kt] = pvl_disc(GHI,Z, doy, varargin)
% PVL_DISC Estimate Direct Normal Irradiance from Global Horizontal Irradiance using the DISC model
%
% Syntax
% [DNI, Kt] = pvl_disc(GHI,Z, doy)
% [DNI, Kt] = pvl_disc(GHI,Z, doy, pressure)

% Description
%   The DISC algorithm converts global horizontal irradiance to direct
%   normal irradiance through empirical relationships between the global
%   and direct clearness indices. 
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
%   pressure - a scalar or vector of values providing the site pressure in
%     Pascal. If pressure is a vector it must be of the same size as all
%     other vector inputs. pressure must be >=0. Pressure may be measured
%     or an average pressure may be calculated from site altitude. If
%     pressure is omitted, standard pressure (101325 Pa) will be used, this
%     is acceptable if the site is near sea level. If the site is not near
%     sea-level, inclusion of a measured or average pressure is highly
%     recommended.
%
% Output:   
%   DNI - the modeled direct normal irradiance in W/m^2 provided by the
%     Direct Insolation Simulation Code (DISC) model. 
%   Kt - Ratio of global to extraterrestrial irradiance on a horizontal
%     plane.
%
% Sources
%
% [1] Maxwell, E. L., "A Quasi-Physical Model for Converting Hourly 
%   Global Horizontal to Direct Normal Insolation", Technical 
%   Report No. SERI/TR-215-3087, Golden, CO: Solar Energy Research 
%   Institute, 1987.
%
% [2] J.W. "Fourier series representation of the position of the sun". 
%   Found at:
%   http://www.mail-archive.com/sundial@uni-koeln.de/msg01050.html on
%   January 12, 2012
%
% See also PVL_DATE2DOY PVL_EPHEMERIS PVL_ALT2PRES PVL_DIRINT PVL_LOUCHE
% PVL_ORGILL_HOLLANDS PVL_REINDL_1 PVL_REINDL_2 PVL_ERBS

p = inputParser;
p.addRequired('GHI', @(x) all(isnumeric(x) & isvector(x) & x>=0));
p.addRequired('Z', @(x) (all(isnumeric(x) & x<=180 & x>=0 & isvector(x))));
p.addRequired('doy', @(x) (all(isnumeric(x) & isvector(x) & x>=1 & x<367)));
p.addOptional('pressure', 101325, @(x) all(isnumeric(x) & x>=0 & isvector(x)));
p.parse(GHI,Z,doy, varargin{:});

pressure = p.Results.pressure;

% It is unclear in the SERI paper describing the DISC model whether Maxwell 
% used the apparent or true zenith angle for calculating airmass. 
% For calculation of the extraterrestrial horizontal radiation
% (later in the algorithm) it seems as if the true zenith angle should be
% used, but I believe that Kasten meant to use an apparent zenith angle in
% his 1965 (or 1966) paper on the determination of air mass which Maxwell
% uses in his algorithm.

% Initialize variables
GHI = GHI.*ones(max([numel(GHI) numel(Z) numel(doy) numel(pressure)]),1);
Z=Z(:);
doy=doy(:);
pressure = pressure(:);
A = zeros(length(GHI),1);
B = zeros(length(GHI),1);
C = zeros(length(GHI),1);

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

DayAngle = 2.*pi.*(doy-1)./365;
re = 1.00011 + 0.034221 .* cos(DayAngle) + (0.00128) .* sin(DayAngle)...
     +0.000719.*cos(2.*DayAngle) + (7.7E-5).*sin(2.*DayAngle);
I0 = re.*1370;
I0h= I0.*cosd(Z);

Ztemp = Z;
Ztemp(Z>87)=87; % Make sure you don't take the root of a negative.
AM = 1./(cosd(Ztemp) + 0.15.*((93.885-Ztemp).^(-1.253))) .* pressure./ 101325;

Kt = GHI./(I0h); % This Z needs to be the true Zenith angle, not apparent (to get extraterrestrial horizontal radiation)
Kt(Kt<0) = 0;

% For Kt > 0.6, set the components of A, B, and C
A(Kt > 0.6) = -5.743 + 21.77.*Kt(Kt>0.6) -27.49.*Kt(Kt>0.6).^2 +11.56.*Kt(Kt>0.6).^3;
B(Kt > 0.6) = 41.4 - 118.5.*Kt(Kt>0.6) + 66.05.*Kt(Kt>0.6).^2 +31.9.*Kt(Kt>0.6).^3;
C(Kt > 0.6) = -47.01 + 184.2.*Kt(Kt>0.6) - 222.*Kt(Kt>0.6).^2 +73.81.*Kt(Kt>0.6).^3;

% For Kt <= 0.6, set the components of A, B, and C
A(Kt <= 0.6) = 0.512 - 1.56.*Kt(Kt <= 0.6) +2.286.*Kt(Kt <= 0.6).^2 - 2.222.*Kt(Kt <= 0.6).^3;
B(Kt <= 0.6) = 0.37 + 0.962.*Kt(Kt <= 0.6);
C(Kt <= 0.6) = -0.28 + 0.932.*Kt(Kt <= 0.6) - 2.048.*Kt(Kt <= 0.6).^2;

delKn = A + B.*exp(C.*AM);

% In some electronic versions of the DISC model paper, the first term
% "0.886" is crossed out and "0.866" is written underneath. However, in the
% DISC model spreadsheet http://rredc.nrel.gov/solar/models/DISC/ on Feb. 1
% 2012, 0.886 is the correct term.
% However, in FORTRAN code obtained from SRRL on Feb. 9 2012, the value of
% 0.866 is used.
% It's 2 vs. 1, and we are using "0.866" since it has 2 "votes".
%
Knc =  0.866 - 0.122.*AM + 0.0121.*AM.^2 - 0.000653.*AM.^3 + 0.000014.*AM.^4;
Kn = Knc - delKn;

%DNI(Kt <=0) = 0;
% DNI(Kt >0 ) = ((Knc(Kt>0)) - delKn(Kt>0)).*Ioh(Kt>0);
% DNI(Z>87 | DNI<0) = 0;

DNI = (Kn).*I0;
DNI(Z>87 | GHI<1 | DNI<0) = 0; % 87 degrees chosen to simulate the FORTRAN code in use by SRRL (from Perez)


% The code above was vectorized for speed based upon the code given in the
% comments below.

% 
% A = zeros(length(GHI),1);
% B = zeros(length(GHI),1);
% C = zeros(length(GHI),1);
% delKn = zeros(length(GHI),1);
% Knc = zeros(length(GHI),1);
% DNI = zeros(length(GHI),1);
% 
% Kt = GHI./(cosd(Z).*Ea);
% 
% for i = 1:length(GHI)
%     if Kt(i) >0.6
%         A(i) = -5.743 + 21.77*Kt(i) -27.49*Kt(i)^2 +11.56*Kt(i)^3;
%         B(i) = 41.4 - 118.5*Kt(i) + 66.06*Kt(i)^2 +31.9*Kt(i)^3;
%         C(i) = -47.01 + 184.2*Kt(i) - 222*Kt(i)^2 +73.81*Kt(i)^3;
%     end
% 
%     if Kt(i) <= 0.6
%         A(i) = 0.512 - 1.56*Kt(i) +2.286*Kt(i)^2 - 2.222*Kt(i)^3;
%         B(i) = 0.37 + 0.962*Kt(i);
%         C(i) = -0.28 + 0.932*Kt(i) - 2.048*Kt(i)^2;
%     end
%     delKn(i) = A(i) + B(i)*exp(C(i).*AM(i));
%     Knc(i) = 0.886 - 0.122*AM(i) + 0.0121*AM(i)^2 - 0.000653*AM(i)^3 + 0.000014*AM(i)^4;
%     if Kt(i) <= 0
%         DNI(i,1) = 0;
%     end
%     if Kt(i) > 0
%         DNI(i) = (Knc(i) - delKn(i))* Ea(i);
%     end
% end
% filter = Z>87 | DNI<0;
% DNI(filter) = 0;
% end