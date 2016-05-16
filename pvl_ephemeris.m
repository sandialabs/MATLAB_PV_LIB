function [SunAz, SunEl, ApparentSunEl, SolarTime]= pvl_ephemeris(Time, Location, varargin)
% PVL_EPHEMERIS Calculates the position of the sun given time, location, and optionally pressure and temperature
%
% Syntax
%   [SunAz, SunEl, ApparentSunEl, SolarTime]=pvl_ephemeris(Time, Location)
%   [SunAz, SunEl, ApparentSunEl, SolarTime]=pvl_ephemeris(Time, Location, Pressure)
%   [SunAz, SunEl, ApparentSunEl, SolarTime]=pvl_ephemeris(Time, Location, Pressure, Temperature)
%   [SunAz, SunEl, ApparentSunEl, SolarTime]=pvl_ephemeris(Time, Location, 'temperature', Temperature)
%
% Description  
%  [SunAz, SunEl, ApparentSunEl, SolarTime]=pvl_ephemeris(Time, Location)
%      Uses the given time and location structs to give sun positions with
%      the pressure assumed to be 1 atm (101325 Pa) and the temperature
%      assumed to be 12 C.
%   [SunAz, SunEl, ApparentSunEl, SolarTime]=pvl_ephemeris(Time, Location, Pressure)
%      Uses the given time and location structs with the given pressure to
%      determine sun positions. The temperature is assumed to be 12C.
%      Pressure must be given in Pascals (1atm = 101325 Pa). If site pressure
%      is unknown but site altitude is known, a conversion function may be
%      used.
%   [SunAz, SunEl, ApparentSunEl, SolarTime]=pvl_ephemeris(Time, Location, Pressure, Temperature)
%      Uses the given time and location structs with the given pressure and
%      temperature to determine sun positions. Pressure must be given in
%      Pascals, and temperature must be given in C.
%   [SunAz, SunEl, ApparentSunEl, SolarTime]=pvl_ephemeris(Time, Location, 'temperature', Temperature)
%      Uses the given time and location structs with the given temperature
%      (in C) to determine sun positions. Default pressure is 101325 Pa.
%
% Input Parameters:
%   Time is a struct with the following elements, note that all elements
%     can be column vectors, but they must all be the same length. Time is
%     entered as a struct which must include a value for the offset from
%     UTC. Time is absolutely specified by the date, time, and the number
%     of hours of offset from that date and time to UTC. For example, if
%     you live in Boston, USA, and your data is timestamped in local standard
%     time, your UTC offset should be -5.
%
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
%   Location is a struct with the following elements 
%   Location.latitude = vector or scalar latitude in decimal degrees (positive is
%     northern hemisphere)
%   Location.longitude = vector or scalar longitude in decimal degrees (positive is 
%     east of prime meridian)
%   Location.altitude = an optional component of the Location struct, not
%     used in the ephemeris code directly, but it may be used to calculate
%     standard site pressure (see pvl_alt2pres function)
% 
% Output Parameters:
%   SunAz = Azimuth of the sun in decimal degrees from North. 0 = North to 270 = West
%   SunEl = Actual elevation (not accounting for refraction)of the sun 
%     in decimal degrees, 0 = on horizon. The complement of the True Zenith
%     Angle.
%   ApparentSunEl = Apparent sun elevation accounting for atmospheric 
%     refraction. This is the complement of the Apparent Zenith Angle.
%   SolarTime = solar time in decimal hours (solar noon is 12.00).
%
% References
%   Grover Hughes' class and related class materials on Engineering 
%   Astronomy at Sandia National Laboratories, 1985.
%
% See also PVL_MAKETIMESTRUCT PVL_MAKELOCATIONSTRUCT PVL_ALT2PRES
%          PVL_GETAOI PVL_SPA

p = inputParser;
p.addRequired('Time',@isstruct);
p.addRequired('Location',@isstruct);
p.addOptional('pressure',101325, @(x) all(isvector(x) & isnumeric(x) & x>=0));
p.addOptional('temperature',12, @(x) all(isvector(x) & isnumeric(x) & x>=-273.15));
p.parse(Time, Location, varargin{:});



Latitude = p.Results.Location.latitude(:);
% the inversion of longitude is due to the fact that this code was
% originally written for the convention that positive longitude were for
% locations west of the prime meridian. However, the correct convention (as
% of 2009) is to use negative longitudes for locations west of the prime
% meridian. Therefore, the user should input longitude values under the
% correct convention (e.g. Albuquerque is at -106 longitude), but it needs
% to be inverted for use in the code.
Longitude = -1 * p.Results.Location.longitude(:);

Year = p.Results.Time.year(:);
Month = p.Results.Time.month(:);
Day = p.Results.Time.day(:);
Hour = p.Results.Time.hour(:);
Minute = p.Results.Time.minute(:);
Second = p.Results.Time.second(:);
% the inversion of UTC offset is due to the fact that this code was
% originally written for the convention that positive offset values were
% positive for locations west of the prime merdian. However, the correct
% convention (as of 2009) is to use negative offset codes for locaitons
% west of the prime merdian. Therefore, the user should input offset values
% under the correct convention (e.g. EST = -5), but it needs to be inverted
% for use in the following code.
TZone = -1 * p.Results.Time.UTCOffset(:); 

if isscalar(p.Results.pressure)
    pressure =  p.Results.pressure*ones(size(Year));
else
    pressure = p.Results.pressure(:);
end

if isscalar(p.Results.temperature)
    temperature = p.Results.temperature*ones(size(Year));
else
    temperature = p.Results.temperature(:);
end

if ~(numel(pressure)==numel(temperature)) && (numel(pressure)==numel(Time.Year)) && ((numel(pressure)==numel(Latitude)) || isscalar(Latitude))
    error(['An error has occurred in ',mfilename,'. One of the input'...
        ' variables is not of correct size.']);
end


DayOfYear = pvl_date2doy(Year, Month, Day);


DecHours = Hour + Minute./60 + Second./3600;

RadtoDeg = 180 / pi;
DegtoRad = pi / 180;

Abber = 20/3600;
LatR = Latitude * DegtoRad;
UnivDate = DayOfYear + floor((DecHours + TZone)/24);
UnivHr = mod((DecHours + TZone), 24);
Yr = Year-1900;
YrBegin = 365 * Yr + floor((Yr-1)/4)-0.5;
Ezero = YrBegin + UnivDate;
T = Ezero / 36525;
GMST0 = 6/24 +38/1440 + (45.836 + 8640184.542 * T + 0.0929 * T.^2)/86400;
GMST0 = 360 * (GMST0 - floor(GMST0));
GMSTi = mod(GMST0 + 360*(1.0027379093 * UnivHr / 24),360);
LocAST = mod((360 + GMSTi - Longitude), 360);
EpochDate = Ezero + UnivHr / 24;
T1 = EpochDate / 36525;
ObliquityR = DegtoRad * (23.452294 - 0.0130125 * T1 - 0.00000164 * T1.^2 ...
    + 0.000000503 * T1.^3);
MlPerigee = 281.22083 + 0.0000470684 * EpochDate + 0.000453 * T1 .^ 2 + ...
    0.000003 * T1 .^ 3;
MeanAnom = mod((358.47583 + 0.985600267 * EpochDate - 0.00015 * T1 .^ 2 - ... 
    0.000003 * T1 .^ 3), 360);
Eccen = 0.01675104 - 0.0000418 * T1 - 0.000000126 * T1 .^ 2;
EccenAnom = MeanAnom;
E=0;

while max(abs(EccenAnom - E)) > 0.0001;
    E = EccenAnom;
    EccenAnom = MeanAnom + RadtoDeg .* Eccen .* sin(DegtoRad .* E);
end

TrueAnom = 2 * mod(RadtoDeg * atan2(((1 + Eccen) ./ (1 - Eccen)).^ 0.5 .* tan(DegtoRad * EccenAnom / 2), 1), 360) ;
EcLon = mod(MlPerigee + TrueAnom, 360) - Abber ;
EcLonR = DegtoRad * EcLon;
DecR = asin(sin(ObliquityR) .* sin(EcLonR));
Dec = RadtoDeg * DecR;
RtAscen = RadtoDeg * atan2(cos(ObliquityR).*(sin(EcLonR)),cos(EcLonR));
HrAngle = LocAST - RtAscen ;
HrAngleR = DegtoRad .* HrAngle ; 

HrAngle = HrAngle - (360 .* sign(HrAngle) .* (abs(HrAngle) > 180));


SunAz = RadtoDeg .* atan2(-1 * sin(HrAngleR), cos(LatR) .* tan(DecR) - sin(LatR) .* cos(HrAngleR));
SunAz = SunAz + (SunAz < 0) * 360; %shift from range of [-180,180] to [0,360]
SunEl = asind(cos(LatR) .* cos(DecR) .* cos(HrAngleR) + sin(LatR) .* sin(DecR));

SolarTime = (180 + HrAngle) / 15;

% Calculate the refraction of the sun until the actual center of the sun is
% 1 degree below the horizon.
TanEl = tan(DegtoRad * SunEl);
Refract = zeros(length(SunEl),1) + ...
    (and(SunEl > 5, SunEl <= 85) .* (58.1 ./ TanEl - 0.07 ./ (TanEl.^3) + .000086 ./ (TanEl.^5))) + ...
    (and(SunEl > -0.575, SunEl <=5) .* (SunEl .* (-518.2 + SunEl .* (103.4 + SunEl .* (-12.79 + SunEl .* 0.711))) +1735)) + ...
    (and(SunEl > -1 ,SunEl <= -0.575) .* ((-20.774 ./ TanEl)));


Refract = Refract .* (283 ./ (273 + temperature)) .* pressure ./ 101325 ./ 3600;

% Generate apparent sun elevation including refraction
ApparentSunEl = SunEl + Refract;




    

    