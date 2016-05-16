function DOY=pvl_date2doy(Year, Month, Day)
% PVL_DATE2DOY Determine day of year from year, month of year, and day of month
%
% Syntax
%   DOY = pvl_date2doy(Year, Month, Day)
% 
% Description
%   Calculate the day of the year given the year, month of year, and day of
%   month in the Gregorian calendar. 
%
% Inputs: 
%       Year - scalar or vector of years (e.g. 1990). Year(i) corresponding 
%           to Month(i) and Day(i). 
%       Month - scalar or vector of month numbers (1-12). Month must be >=1 and <13
%           Month must be of the same length as Year.
%       Day - scalar or vector of day of month numbers. Must be of same
%           length as Year. Day must be >=1 and <32. Day may either be a 
%           whole or fractional day (e.g. the 12th of the month at noon 
%           could be 12.5); if fractional, pvl_date2doy will maintain the fractional day.
%
%   If Year or Month are non-integers, the decimal portion of the values 
%   will be dropped (floor function). The calculation utilizes the 400 year
%   cycle for leap years.
%
% Outputs:
%   DOY - A column vector of the day of year (whole or
%   fractional) as a floating point number from 1 to 366.999... If any of
%   the input values are NaN, the corresponding output value will be NaN.
%
%   Note that while Day must be >=1 and <32. The calculation does NOT check
%   to ensure that the day of month is valid for a given month. Thus
%   pvl_date2doy(2012, 4, 31) = pvl_date2doy(2012, 5, 1). 
%
% See also PVL_LEAPYEAR PVL_DOY2DATE DATEVEC

p=inputParser;
p.addRequired('Year', @(x) all((isnumeric(x) | isnan(x)) & isvector(x)));
p.addRequired('Month', @(x) all(isvector(x) & (isnan(x) | (isnumeric(x) & x>=1 & x<13 & (numel(x)==numel(Year))))));
p.addRequired('Day', @(x) all(isvector(x) & (isnan(x) | (isnumeric(x) & x>=1 & x<32 & (numel(x)==numel(Year))))));
p.parse(Year, Month, Day);

Month = floor(Month(:));
Year = floor(Year(:));
Day = Day(:);

% This error check now taken care of with the input parser
% if (length(Year) ~= length(Month)) || (length(Year) ~= length(Day))
%     error('Error in dayofyear. Year, Month, and Day must be of same length.')
% end

% Determine if Year is a leap year.
LY = pvl_leapyear(Year);

% Calculate day of year for perpetual (non-leap) years
MonthOffset = [1;2;0;1;1;2;2;3;4;4;5;5]; %Offset from (n-1)*30
Month(isnan(Month))=1; % Can't use NaN as an index, so set all NaNs to 1. Take care of the NaNs in output later.
MonthStart = ((Month - 1) .* 30) + MonthOffset(Month);

% Day of Year is the day of year of the first of the month in a perpetual
% year, plus the day of month less 1, plus 1 if the month is greater 
% than 2 AND it is a leapyear
DOY = MonthStart + (Day - 1) + (Month > 2).*LY;

% If any of the input were NaN, make all corresponding outputs NaN
DOY(isnan(p.Results.Year) | isnan(p.Results.Month) | isnan(p.Results.Day)) = NaN;