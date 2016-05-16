function [SunAz, SunEl, ApparentSunEl]= pvl_spa(Time, Location, varargin)
% PVL_SPA Calculates the position of the sun given time, location, and optionally pressure and temperature
%
% Syntax
%   [SunAz, SunEl, ApparentSunEl]=pvl_spa(Time, Location)
%   [SunAz, SunEl, ApparentSunEl]=pvl_spa(Time, Location, Pressure)
%   [SunAz, SunEl, ApparentSunEl]=pvl_spa(Time, Location, Pressure, Temperature)
%   [SunAz, SunEl, ApparentSunEl]=pvl_spa(Time, Location, Pressure, Temperature, delta_t)
%   [SunAz, SunEl, ApparentSunEl]=pvl_spa(Time, Location, 'temperature', Temperature)
%   [SunAz, SunEl, ApparentSunEl]=pvl_spa(Time, Location, 'delta_t', delta_t)
%
% Description  
%  [SunAz, SunEl, ApparentSunEl]=pvl_spa(Time, Location)
%      Uses the given time and location structs to give sun positions with
%      the pressure assumed to be 1 atm (101325 Pa) and the temperature
%      assumed to be 12 C.
%   [SunAz, SunEl, ApparentSunEl]=pvl_spa(Time, Location, Pressure)
%      Uses the given time and location structs with the given pressure to
%      determine sun positions. The temperature is assumed to be 12C.
%      Pressure must be given in Pascals (1atm = 101325 Pa). If site pressure
%      is unknown but site altitude is known, a conversion function may be
%      used.
%   [SunAz, SunEl, ApparentSunEl]=pvl_spa(Time, Location, Pressure, Temperature)
%      Uses the given time and location structs with the given pressure and
%      temperature to determine sun positions. Pressure must be given in
%      Pascals, and temperature must be given in C.
%   [SunAz, SunEl, ApparentSunEl]=pvl_spa(Time, Location, 'temperature', Temperature)
%      Values 'temperature', 'pressure', and 'delta_t', may also be entered
%      in any configuration using the Parameter/Value method of variable
%      entry. Simply use the quoted name of the parameter, then a comma,
%      then the parameter's value.
%
% Input Parameters:
%   Time is a struct with the following elements, note that all elements
%     can be column vectors, but they must all be the same length
%   Time.year = The year in the gregorian calendar
%   Time.month = the month of the year (January = 1 to December = 12)
%   Time.day = the day of the month
%   Time.hour = the hour of the day
%   Time.minute = the minute of the hour
%   Time.second = the second of the minute
%   Time.UTCOffset = the UTC offset code, using the convention
%     that a positive UTC offset is for time zones east of the prime meridian
%     (e.g. EST = -5)
% 
%   Location is a struct with the following elements:
%   Location.latitude = vector or scalar latitude in decimal degrees (positive is
%     northern hemisphere)
%   Location.longitude = vector or scalar longitude in decimal degrees (positive is 
%     east of prime meridian)
%   Location.altitude = vector or scalar site elevation (or altitude) in
%     meters. Required for the SPA calculation.
%
%   Pressure = site pressure in Pascal. This value is optional, but if it
%      is omitted, the default value is 101325 (sea level standard
%      pressure) is used.
%
%   Temperature = site temperature in C. This input is optional, and if
%     omitted, a value of 12C is used.
%
%   delta_t = The current difference, in seconds, between international atomic time
%      (TAI) and UT1. This value is obtained by observation (see [2]). If
%      omitted, the default value is 66.3 + 0.6175*(Year-2012). This
%      default prediction can be changed as necessary.
% 
% Output Parameters:
%   SunAz = Azimuth of the sun in decimal degrees from North. 0 = North to 270 = West
%   SunEl = Actual elevation (not accounting for refraction)of the sun 
%     in decimal degrees, 0 = on horizon. The complement of the True Zenith
%     Angle.
%   ApparentSunEl = Apparent sun elevation accounting for atmospheric 
%     refraction. This is the complement of the Apparent Zenith Angle.
%
% References
%   [1] Reda, I. and Andreas, A. Solar Position Algorithm for Solar
%       Radiaiton Applications. NREL Report No. TP-560-34302. Revised
%       January 2008. http://www.nrel.gov/docs/fy08osti/34302.pdf as viewed
%       2012-05-07
%   [2] DeltaT projections and current values from US Naval Observatory,
%       http://asa.usno.navy.mil/SecK/DeltaT.html as viewed 2012-05-07
%
% See also PVL_MAKETIMESTRUCT PVL_MAKELOCATIONSTRUCT PVL_ALT2PRES
%          PVL_GETAOI PVL_EPHEMERIS

p = inputParser;
p.addRequired('Time',@isstruct);
p.addRequired('Location',@isstruct);
p.addOptional('pressure',101325, @(x) all(isvector(x) & isnumeric(x) & x>=0));
p.addOptional('temperature',12, @(x) all(isvector(x) & isnumeric(x) & x>=-273.15));
p.addOptional('delta_t', 66.3 + .6175*(Time.year-2012), @(x) all(isvector(x) & isnumeric(x) & x>-8000 & x<8000));
p.parse(Time, Location, varargin{:});

% Read in optional values from the struct p
delta_t = p.Results.delta_t;
pressure = p.Results.pressure ./ 100; % Convert to millibars from Pa
temperature = p.Results.temperature;

% For application of atmospheric refraction correction
SunRadius = 0.26667; % Apparent Sun Radius in degrees
AtmosRefract = 0.5667; % Refraction of the sun on the horizon in degrees


% Calculate the year, month, and day. If month <=2, then reduce the year by
% 1 and increase the month by 12.
Y = Time.year - 1 .* (Time.month <= 2);
M = Time.month + 12 .* (Time.month <=2);
D = Time.day + (Time.hour - Time.UTCOffset)./24 + Time.minute ./ (24*60) + Time.second ./ (24 * 60 *60);

% Note: Here I deviate from the SPA algorithm as defined in the 2008 SPA
% paper. In the paper they specify using an "INT" function such as int64.
% However, in the C code, they use a "floor" function instead of
% typecasting. I believe that the "floor" function is correct, since using
% a typecasting approach can give identical JD values for different days
% (e.g. Jan 1, -122 and Dec, 31 -123).
A = fix(Y./100);
B = (2 - A + fix(A/4));


% Calculate Julian Day and other time metrics
JD = (fix(365.25 .* (Y + 4716))) + (fix(30.6001 .* (M + 1))) + D + 0 - 1524.5;

% There may be errors in calculating leap years for Julian calendar dates
% before year 1. Table A4.1 in the SPA code lists Feb. 29, -1000 as a valid
% day, but MATLAB's datenum function condsiders Feb. 29, -1000 and Mar. 1, -1000 
% as the same day (which means it believes that -1000 is NOT a leap year).
% The equations here should work correctly given the Feb. 29, -1000 date,
% as long as the MATLAB function datenum is not used previously.
if any(JD < 1721421.5)
    warning(['Dates prior to Julian year 1 detected in ', mfilename, ...
        '. Please note that leap-years before year 1 may be calculated',...
        ' incorrectly.']);
end

% Add B to JD only if JD>2299160
JD(JD>2299160)=JD(JD>2299160)+B(JD>2299160);
JDE = JD + delta_t/86400;
JC = (JD - 2451545) / 36525;
JCE = (JDE - 2451545) / 36525;
JME = JCE/10;

JME = JME(:); % Ensure column vector

% Attempt to clear a little bit of memory
clear Time;
%% 
% Get the required tables (Appendix A.4 in SPA paper)
[L0table, L1table, L2table, L3table, L4table, L5table,...
    B0table, B1table,...
    R0table, R1table, R2table, R3table, R4table,...
    Ytable, Psitable, Epsilontable...
    ] = get_tables();
%% Calculate the Earth heliocentric longitude 3.2.1 - 3.2.6
Asub0 = L0table(:,1)';
Bsub0 = L0table(:,2)';
Csub0 = L0table(:,3)';
L0 = sum((ones(size(JME)) * Asub0) .* cos((ones(size(JME)) * Bsub0) + (JME * Csub0)), 2);

Asub1 = L1table(:,1)';
Bsub1 = L1table(:,2)';
Csub1 = L1table(:,3)';
L1 = sum((ones(size(JME)) * Asub1) .* cos((ones(size(JME)) * Bsub1) + (JME * Csub1)), 2);

Asub2 = L2table(:,1)';
Bsub2 = L2table(:,2)';
Csub2 = L2table(:,3)';
L2 = sum((ones(size(JME)) * Asub2) .* cos((ones(size(JME)) * Bsub2) + (JME * Csub2)), 2);

Asub3 = L3table(:,1)';
Bsub3 = L3table(:,2)';
Csub3 = L3table(:,3)';
L3 = sum((ones(size(JME)) * Asub3) .* cos((ones(size(JME)) * Bsub3) + (JME * Csub3)), 2);

Asub4 = L4table(:,1)';
Bsub4 = L4table(:,2)';
Csub4 = L4table(:,3)';
L4 = sum((ones(size(JME)) * Asub4) .* cos((ones(size(JME)) * Bsub4) + (JME * Csub4)), 2);

Asub5 = L5table(:,1);
Bsub5 = L5table(:,2);
Csub5 = L5table(:,3);
L5 = sum((ones(size(JME)) * Asub5) .* cos((ones(size(JME)) * Bsub5) + (JME * Csub5)), 2);

L_rad = 1E-8 * (L0 + L1.*JME + L2.*(JME.^2) +L3.*(JME.^3) + L4.*(JME.^4) + L5.*(JME.^5));
L_deg = L_rad * 180 / pi;
F_L = L_deg/360 - floor(L_deg/360);
L_deg(L_deg < 0) = 360 - 360 .* F_L(L_deg < 0);
L_deg(L_deg >= 0) = 360 .* F_L(L_deg >= 0);

% Attempt to clear a little bit of memory
clear L_rad;
clear F_L;
clear L0; clear L1; clear L2; clear L3; clear L4; clear L5;
%% Calculate the Earth heliocentric latitude 3.2.7
Asub0 = B0table(:,1)';
Bsub0 = B0table(:,2)';
Csub0 = B0table(:,3)';
B0 = sum((ones(size(JME)) * Asub0) .* cos((ones(size(JME)) * Bsub0) + (JME * Csub0)), 2);

Asub1 = B1table(:,1)';
Bsub1 = B1table(:,2)';
Csub1 = B1table(:,3)';
B1 = sum((ones(size(JME)) * Asub1) .* cos((ones(size(JME)) * Bsub1) + (JME * Csub1)), 2);

B2 = 0;

B3 = 0;

B4 = 0;

B5 = 0;

B_rad = 1E-8 * (B0 + B1.*JME + B2.*(JME.^2) +B3.*(JME.^3) + B4.*(JME.^4) + B5.*(JME.^5));
B_deg = B_rad * 180 / pi;
% 
% F_B = B_deg/360 - floor(B_deg/360);
% B_deg(B_deg < 0) = 360 - 360 .* F_B(B_deg < 0);
% B_deg(B_deg >= 0) = 360 .* F_B(B_deg >= 0);

% Attempt to clear a little bit of memory
clear B_rad;
clear B0; clear B1;

%% Calculate the Earth radius vector 3.2.8
Asub0 = R0table(:,1)';
Bsub0 = R0table(:,2)';
Csub0 = R0table(:,3)';
R0 = sum((ones(size(JME)) * Asub0) .* cos((ones(size(JME)) * Bsub0) + (JME * Csub0)), 2);

Asub1 = R1table(:,1)';
Bsub1 = R1table(:,2)';
Csub1 = R1table(:,3)';
R1 = sum((ones(size(JME)) * Asub1) .* cos((ones(size(JME)) * Bsub1) + (JME * Csub1)), 2);

Asub2 = R2table(:,1)';
Bsub2 = R2table(:,2)';
Csub2 = R2table(:,3)';
R2 = sum((ones(size(JME)) * Asub2) .* cos((ones(size(JME)) * Bsub2) + (JME * Csub2)), 2);

Asub3 = R3table(:,1)';
Bsub3 = R3table(:,2)';
Csub3 = R3table(:,3)';
R3 = sum((ones(size(JME)) * Asub3) .* cos((ones(size(JME)) * Bsub3) + (JME * Csub3)), 2);

Asub4 = R4table(:,1)';
Bsub4 = R4table(:,2)';
Csub4 = R4table(:,3)';
R4 = sum((ones(size(JME)) * Asub4) .* cos((ones(size(JME)) * Bsub4) + (JME * Csub4)), 2);

R5 = 0;

R = 1E-8 * (R0 + R1.*JME + R2.*(JME.^2) +R3.*(JME.^3) + R4.*(JME.^4) + R5.*(JME.^5));

% Attempt to clear a little bit of memory
clear R0; clear R1; clear R2; clear R3; clear R4; 

%% Calculate the geocentric longitude and latitude 3.3
Theta_deg = L_deg + 180;
F_Theta = Theta_deg/360 - floor(Theta_deg/360);
Theta_deg(Theta_deg < 0) = 360 - 360 .* F_Theta(Theta_deg < 0);
Theta_deg(Theta_deg >= 0) = 360 .* F_Theta(Theta_deg >= 0);

Beta_deg = -1.* B_deg;

% Attempt to clear a little bit of memory
clear F_Theta;

%% Calculate the nutation in longitude and obliquity 3.4
X0 = 297.85036 + 445267.111480 .* JCE - .0019142 .* JCE.^2 + 1/189474 .* JCE.^3;
X1 = 357.52772 + 35999.050340 .* JCE - .0001603 .* JCE.^2 - 1/300000 .* JCE.^3;
X2 = 134.96298 + 477198.867398 .* JCE + .0086972 .* JCE.^2 + 1/56250 .* JCE.^3;
X3 = 93.27191 + 483202.017538 .* JCE - .0036825 .* JCE.^2 + 1/327270 .* JCE.^3;
X4 = 125.04452 - 1934.136261 .* JCE + .0020708 .* JCE.^2 + 1/450000 .* JCE.^3;
X = [X0 , X1 , X2 , X3 , X4];

sumXY = Ytable * X';

deltaPsi_i = ((Psitable(:,1) * ones(1,length(JCE))) + ... % (a +
    ((Psitable(:,2) * ones(1,length(JCE))) .* (ones(length(Psitable(:,2)),1) * JCE'))) .* ... % b * JCE') * 
    sind(sumXY); % sind(Sum(Xj * Yi,j))
    
deltaEpsilon_i = ((Epsilontable(:,1) * ones(1,length(JCE))) + ... % (c +
    ((Epsilontable(:,2) * ones(1,length(JCE))) .* (ones(length(Epsilontable(:,2)),1) * JCE'))) .* ... % d * JCE') * 
    cosd(sumXY); % cosd(Sum(Xj * Yi,j))

deltaPsi_deg = (sum(deltaPsi_i, 1)/36000000)';
deltaEpsilon_deg = (sum(deltaEpsilon_i, 1)/36000000)';

% Attempt to clear a little bit of memory
clear deltaPsi_i;
clear deltaEpsilon_i;
clear sumXY;
clear X; clear X0; clear X1; clear X2; clear X3; clear X4;
%% Calculate the true obliquity of the ecliptic 3.5

U = JME ./ 10;

% Note, in the 2008 version of the SPA paper and C code, there is a
% discrepency between the stated algorithm and the the C code for this
% equation. The algorithm (eqn 24) has the second term as -4680.93, while
% the C code given in the appendix uses -4680.96. However, the actual .c
% code which can be obtained via downlaod uses -4680.93, so that is used
% here.
Epsilon0_deg = 84381.448 - 4680.93 .* U - 1.55 .* U.^2 + 1999.25 .* U.^3 - ...
    51.38 .* U.^4 - 249.67 .* U.^5 - 39.05.* U.^6 + 7.125 .* U.^7 + ...
    27.87 .* U.^8 + 5.79 .* U.^9 + 2.45 .* U.^10;

Epsilon_deg = Epsilon0_deg./3600 + deltaEpsilon_deg;

% Attempt to clear a little bit of memory
clear U;
%% Calculate the aberration correction 3.6
deltaTau_deg = -20.4898./ (3600 .* R);

%% Calculate the apparent sun longitude 3.7
Lambda_deg = Theta_deg + deltaPsi_deg + deltaTau_deg;

%% Calculate the apparent sidereal time at Greenwich at any given time 3.8
nu0_deg = 280.46061837 + 360.98564736629.* (JD - 2451545) + ...
    0.000387933 .* JC.^2 - (JC.^3)/38710000;

F_nu0 = nu0_deg/360 - floor(nu0_deg/360);
nu0_deg(nu0_deg < 0) = 360 - 360 .* F_nu0(nu0_deg < 0);
nu0_deg(nu0_deg >= 0) = 360 .* F_nu0(nu0_deg >= 0);


nu_deg = nu0_deg + deltaPsi_deg .* cosd(Epsilon_deg);

% Attempt to clear a little bit of memory
clear F_nu0;
%% Calculate the geocentric sun right ascension 3.9
Alpha_rad = atan2(sind(Lambda_deg) .* cosd(Epsilon_deg) - ...
    tand(Beta_deg) .* sind(Epsilon_deg),...
    cosd(Lambda_deg));
Alpha_deg = Alpha_rad .* 180 ./ pi;


F_Alpha = Alpha_deg/360 - floor(Alpha_deg/360);
Alpha_deg(Alpha_deg < 0) = 360 - 360 .* F_Alpha(Alpha_deg < 0);
Alpha_deg(Alpha_deg >= 0) = 360 .* F_Alpha(Alpha_deg >= 0);

% Attempt to clear a little bit of memory
clear F_Alpha; clear Alpha_rad;
%% Calculate the geocentric sun declination 3.10

Delta_deg = asind(sind(Beta_deg) .* cosd(Epsilon_deg) + ...
    cosd(Beta_deg) .* sind(Epsilon_deg) .* sind(Lambda_deg));

%% Calculate the observer local hour angle 3.11

H_deg = nu_deg + Location.longitude - Alpha_deg;

F_H = H_deg/360 - floor(H_deg/360);
H_deg(H_deg < 0) = 360 - 360 .* F_H(H_deg < 0);
H_deg(H_deg >= 0) = 360 .* F_H(H_deg >= 0);

% Attempt to clear a little bit of memory
clear F_H;
%% Calculate the topocentric sun right ascension 3.12

Xi_deg = 8.794 ./ (3600 .* R);

u_deg = atand(0.99664719 .* tand(Location.latitude));

x = cosd(u_deg) + Location.altitude./6378140 .* cosd(Location.latitude);

y = 0.99664719 .* sind(u_deg) + Location.altitude ./ 6378140 .*sind(Location.latitude);

deltaAlpha_rad = atan2(-1.*x .* sind(Xi_deg) .* sind(H_deg),...
    cosd(Delta_deg) + (-1.*x .* sind(Xi_deg) .* cosd(H_deg)));

deltaAlpha_deg = deltaAlpha_rad .* 180 ./ pi;

Alphaprime_deg = Alpha_deg + deltaAlpha_deg;

Deltaprime_rad = atan2( (sind(Delta_deg)-1.*y.*sind(Xi_deg).*cosd(deltaAlpha_deg)) , ...
    cosd(Delta_deg)-x.*sind(Xi_deg).*cosd(H_deg));
Deltaprime_deg = Deltaprime_rad .* 180 ./ pi;

% Attempt to clear a little bit of memory
clear Xi_deg; clear u_deg; clear x; clear y;
%% Calculate the topocentric local hour angle 3.13

Hprime_deg = H_deg - deltaAlpha_deg;

%% Calculate the topocentric zenith angle 3.14

e0_deg = asind(sind(Location.latitude).*sind(Deltaprime_deg) + ...
    cosd(Location.latitude).*cosd(Deltaprime_deg) .* cosd(Hprime_deg));


deltae_deg = pressure/1010 .* 283 ./ (273 + temperature) .* ...
    1.02 ./ ((60 .* tand(e0_deg + (10.3./ (e0_deg+5.11)))));

ApplyRefractionFilter = (e0_deg >= -1*(AtmosRefract+SunRadius));

e_deg = e0_deg;
e_deg(ApplyRefractionFilter)=e_deg(ApplyRefractionFilter)+ deltae_deg(ApplyRefractionFilter);

apparentzenith_deg = 90 - e_deg;

SunEl = 90-(90-e0_deg);
ApparentSunEl = 90-(apparentzenith_deg);

%% Calculate the topocentric azimuth angle 3.15

Gamma_rad = atan2(sind(Hprime_deg) , ...
    cosd(Hprime_deg).*sind(Location.latitude) - tand(Deltaprime_deg).*cosd(Location.latitude));

Gamma_deg = Gamma_rad * 180/pi;

F_Gamma = Gamma_deg/360 - floor(Gamma_deg/360);
Gamma_deg(Gamma_deg < 0) = 360 - 360 .* F_Gamma(Gamma_deg < 0);
Gamma_deg(Gamma_deg >= 0) = 360 .* F_Gamma(Gamma_deg >= 0);

azimuth_deg = Gamma_deg + 180;
F_azimuth_deg = azimuth_deg/360 - floor(azimuth_deg/360);
azimuth_deg(azimuth_deg < 0) = 360 - 360 .* F_azimuth_deg(azimuth_deg < 0);
azimuth_deg(azimuth_deg >= 0) = 360 .* F_azimuth_deg(azimuth_deg >= 0);
SunAz = azimuth_deg;

%% A.1
SunMeanLong_deg = 280.4664567 + 360007.6982779 .* JME + 0.03032028 .* JME .^2 ...
    + 1./49931.*JME.^3 - 1./15300.*JME.^4 - 1./2E6 .*JME.^5;

F_SunMeanLong_deg = SunMeanLong_deg/360 - floor(SunMeanLong_deg/360);
SunMeanLong_deg(SunMeanLong_deg < 0) = 360 - 360 .* F_SunMeanLong_deg(SunMeanLong_deg < 0);
SunMeanLong_deg(SunMeanLong_deg >= 0) = 360 .* F_SunMeanLong_deg(SunMeanLong_deg >= 0);

E_deg = SunMeanLong_deg - 0.0057183 - Alpha_deg + deltaPsi_deg .* cosd(Epsilon_deg);

F_E_deg = E_deg/360 - floor(E_deg/360);
E_deg(E_deg < 0) = 360 - 360 .* F_E_deg(E_deg < 0);
E_deg(E_deg >= 0) = 360 .* F_E_deg(E_deg >= 0);

E_minutes = E_deg .* 4;

E_minutes(E_minutes < -20) = E_minutes(E_minutes < -20) + 1440;
E_minutes(E_minutes > 20) = E_minutes(E_minutes > 20) - 1440;

%% A.2


end
%%
% This function just loads up the tables necessary
function [L0table, L1table, L2table, L3table, L4table, L5table,...
    B0table, B1table,...
    R0table, R1table, R2table, R3table, R4table,...
    Ytable, Psitable, Epsilontable...
    ] = get_tables()
L0table = ...
[175347046.0,0,0;
 3341656.0,4.6692568,6283.07585;
 34894.0,4.6261,12566.1517;
 3497.0,2.7441,5753.3849;
 3418.0,2.8289,3.5231;
 3136.0,3.6277,77713.7715;
 2676.0,4.4181,7860.4194;
 2343.0,6.1352,3930.2097;
 1324.0,0.7425,11506.7698;
 1273.0,2.0371,529.691;
 1199.0,1.1096,1577.3435;
 990,5.233,5884.927;
 902,2.045,26.298;
 857,3.508,398.149;
 780,1.179,5223.694;
 753,2.533,5507.553;
 505,4.583,18849.228;
 492,4.205,775.523;
 357,2.92,0.067;
 317,5.849,11790.629;
 284,1.899,796.298;
 271,0.315,10977.079;
 243,0.345,5486.778;
 206,4.806,2544.314;
 205,1.869,5573.143;
 202,2.458,6069.777;
 156,0.833,213.299;
 132,3.411,2942.463;
 126,1.083,20.775;
 115,0.645,0.98;
 103,0.636,4694.003;
 102,0.976,15720.839;
 102,4.267,7.114;
 99,6.21,2146.17;
 98,0.68,155.42;
 86,5.98,161000.69;
 85,1.3,6275.96;
 85,3.67,71430.7;
 80,1.81,17260.15;
 79,3.04,12036.46;
 75,1.76,5088.63;
 74,3.5,3154.69;
 74,4.68,801.82;
 70,0.83,9437.76;
 62,3.98,8827.39;
 61,1.82,7084.9;
 57,2.78,6286.6;
 56,4.39,14143.5;
 56,3.47,6279.55;
 52,0.19,12139.55;
 52,1.33,1748.02;
 51,0.28,5856.48;
 49,0.49,1194.45;
 41,5.37,8429.24;
 41,2.4,19651.05;
 39,6.17,10447.39;
 37,6.04,10213.29;
 37,2.57,1059.38;
 36,1.71,2352.87;
 36,1.78,6812.77;
 33,0.59,17789.85;
 30,0.44,83996.85;
 30,2.74,1349.87;
 25,3.16,4690.48];

L1table = ...
 [628331966747.0,0,0;
 206059.0,2.678235,6283.07585;
 4303.0,2.6351,12566.1517;
 425.0,1.59,3.523;
 119.0,5.796,26.298;
 109.0,2.966,1577.344;
 93,2.59,18849.23;
 72,1.14,529.69;
 68,1.87,398.15;
 67,4.41,5507.55;
 59,2.89,5223.69;
 56,2.17,155.42;
 45,0.4,796.3;
 36,0.47,775.52;
 29,2.65,7.11;
 21,5.34,0.98;
 19,1.85,5486.78;
 19,4.97,213.3;
 17,2.99,6275.96;
 16,0.03,2544.31;
 16,1.43,2146.17;
 15,1.21,10977.08;
 12,2.83,1748.02;
 12,3.26,5088.63;
 12,5.27,1194.45;
 12,2.08,4694;
 11,0.77,553.57;
 10,1.3,6286.6;
 10,4.24,1349.87;
 9,2.7,242.73;
 9,5.64,951.72;
 8,5.3,2352.87;
 6,2.65,9437.76;
 6,4.67,4690.48];
  
L2table = [...
 52919 0 0  
 8720 1.0721 6283.0758  
 309 0.867 12566.152  
 27 0.05 3.52  
 16 5.19 26.3  
 16 3.68 155.42  
 10 0.76 18849.23  
 9 2.06 77713.77  
 7 0.83 775.52  
 5 4.66 1577.34  
 4 1.03 7.11  
 4 3.44 5573.14  
 3 5.14 796.3  
 3 6.05 5507.55  
 3 1.19 242.73  
 3 6.12 529.69  
 3 0.31 398.15  
 3 2.28 553.57  
 2 4.38 5223.69  
 2 3.75 0.98];

L3table =[...
 289 5.844 6283.076  
 35 0 0  
 17 5.49 12566.15  
 3 5.2 155.42  
 1 4.72 3.52  
 1 5.3 18849.23  
 1 5.97 242.73]; 
 
L4table = [...
 114 3.142 0  
 8 4.13 6283.08  
 1 3.84 12566.15];

L5table = [1 3.14 0]; 

B0table = ...
 [280.0,3.199,84334.662;
 102.0,5.422,5507.553;
 80,3.88,5223.69;
 44,3.7,2352.87;
 32,4,1577.34];

  
B1table = [...
 9 3.9 5507.55  
 6 1.73 5223.69]; 

R0table = ...
[100013989.0,0,0;
 1670700.0,3.0984635,6283.07585;
 13956.0,3.05525,12566.1517;
 3084.0,5.1985,77713.7715;
 1628.0,1.1739,5753.3849;
 1576.0,2.8469,7860.4194;
 925.0,5.453,11506.77;
 542.0,4.564,3930.21;
 472.0,3.661,5884.927;
 346.0,0.964,5507.553;
 329.0,5.9,5223.694;
 307.0,0.299,5573.143;
 243.0,4.273,11790.629;
 212.0,5.847,1577.344;
 186.0,5.022,10977.079;
 175.0,3.012,18849.228;
 110.0,5.055,5486.778;
 98,0.89,6069.78;
 86,5.69,15720.84;
 86,1.27,161000.69;
 65,0.27,17260.15;
 63,0.92,529.69;
 57,2.01,83996.85;
 56,5.24,71430.7;
 49,3.25,2544.31;
 47,2.58,775.52;
 45,5.54,9437.76;
 43,6.01,6275.96;
 39,5.36,4694;
 38,2.39,8827.39;
 37,0.83,19651.05;
 37,4.9,12139.55;
 36,1.67,12036.46;
 35,1.84,2942.46;
 33,0.24,7084.9;
 32,0.18,5088.63;
 32,1.78,398.15;
 28,1.21,6286.6;
 28,1.9,6279.55;
 26,4.59,10447.39];
   
R1table = ...
 [103019.0,1.10749,6283.07585;
 1721.0,1.0644,12566.1517;
 702.0,3.142,0;
 32,1.02,18849.23;
 31,2.84,5507.55;
 25,1.32,5223.69;
 18,1.42,1577.34;
 10,5.91,10977.08;
 9,1.42,6275.96;
 9,0.27,5486.78];
  
R2table = ...
[4359.0,5.7846,6283.0758;
 124.0,5.579,12566.152;
 12,3.14,0;
 9,3.63,77713.77;
 6,1.87,5573.14;
 3,5.47,18849.23];
 
R3table = [...
 145 4.273 6283.076  
 7 3.92 12566.15];

R4table = [4 2.56 6283.08];

Ytable =  [...
 0 0 0 0 1  
 -2 0 0 2 2  
 0 0 0 2 2  
 0 0 0 0 2  
 0 1 0 0 0  
 0 0 1 0 0  
 -2 1 0 2 2  
 0 0 0 2 1  
 0 0 1 2 2  
 -2 -1 0 2 2  
 -2 0 1 0 0  
 -2 0 0 2 1  
 0 0 -1 2 2  
 2 0 0 0 0  
 0 0 1 0 1  
 2 0 -1 2 2  
 0 0 -1 0 1  
 0 0 1 2 1  
 -2 0 2 0 0  
 0 0 -2 2 1  
 2 0 0 2 2  
 0 0 2 2 2  
 0 0 2 0 0  
 -2 0 1 2 2  
 0 0 0 2 0  
 -2 0 0 2 0  
 0 0 -1 2 1  
 0 2 0 0 0  
 2 0 -1 0 1  
 -2 2 0 2 2  
 0 1 0 0 1  
 -2 0 1 0 1  
 0 -1 0 0 1  
 0 0 2 -2 0  
 2 0 -1 2 1  
 2 0 1 2 2  
 0 1 0 2 2  
 -2 1 1 0 0  
 0 -1 0 2 2  
 2 0 0 2 1  
 2 0 1 0 0  
 -2 0 2 2 2  
 -2 0 1 2 1  
 2 0 -2 0 1  
 2 0 0 0 1  
 0 -1 1 0 0  
 -2 -1 0 2 1  
 -2 0 0 0 1  
 0 0 2 2 1  
 -2 0 2 0 1  
 -2 1 0 2 1  
 0 0 1 -2 0  
 -1 0 1 0 0  
 -2 1 0 0 0  
 1 0 0 0 0  
 0 0 1 2 0  
 0 0 -2 2 2  
 -1 -1 1 0 0  
 0 1 1 0 0  
 0 -1 1 2 2  
 2 -1 -1 2 2  
 0 0 3 2 2  
 2 -1 0 2 2];

Psitable = [-171996 -174.2
-13187 -1.6
-2274 -0.20
2062 0.20
1426 -3.4
712 0.10
-517 1.2
-386 -0.40
-301 0
217 -0.50
-158 0
129 0.10
123 0
63 0
63 0.10
-59 0
-58 -0.1
-51 0
48 0
46 0
-38 0
-31 0
29 0
29 0
26 0
-22 0
21 0
17 -0.1
16 0
-16 0.1
-15 0
-13 0
-12 0
11 0
-10 0
-8 0
7 0
-7 0
-7 0
-7 0
6 0
6 0
6 0
-6 0
-6 0
5 0
-5 0
-5 0
-5 0
4 0
4 0
4 0
-4 0
-4 0
-4 0
3 0
-3 0
-3 0
-3 0
-3 0
-3 0
-3 0
-3 0];

Epsilontable = [92025 8.9
5736 -3.1
977 -0.50
-895 0.50
54 -0.10
-7 0
224 -0.60
200 0
129 -0.10
-95 0.30
0 0
-70 0
-53 0
0 0
-33 0
26 0
32 0
27 0
0 0
-24 0
16 0
13 0
0 0
-12 0
0 0
0 0
-10 0
0 0
-8 0
7 0
9 0
7 0
6 0
0 0
5 0
3 0
-3 0
0 0
3 0
3 0
0 0
-3 0
-3 0
3 0
3 0
0 0
3 0
3 0
3 0
0 0
0 0
0 0
0 0
0 0
0 0
0 0
0 0
0 0
0 0
0 0
0 0
0 0
0 0];

end