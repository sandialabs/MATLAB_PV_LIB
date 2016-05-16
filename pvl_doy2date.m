function [yr, mo, da] = pvl_doy2date(year,dayofyear)
% PVL_DOY2DATE Determines the Year, Month of year, and Day of month given the year and day of year
%
% Syntax
%   [yr mo da] = pvl_doy2date(year, dayofyear)
%
% Description
%   Produces the Year, Month of year, and Day of month in the Gregorian
%   calendar, when given the year and day of year. 
% 
%   Input "year" must be a numeric vector or scalar of gregorian years. 
%   If "year" is a scalar, it may not be NaN.
%   If "year" is not an integer, the fractional part will be removed (floor).
%   Input "dayofyear" must be a numeric vector with each element >=1 and
%   <366 (< 367 for leap years); "dayofyear" may also be NaN. 
%
%   Output "yr" is a column vector of same size as dayofyear returning the 
%       year input. If the input "year" was a scalar, then "yr" is a vector 
%       filled with "year".
%   Output "mo" is a column vector of month of the year (1 = January, 
%       12 = December)
%   Output "da" is a column vector of day of the month, including any 
%       fractional days given in dayofyear.
%   All output variables are set to NaN if any of the corresponding input
%       values are NaN.
%
% See also PVL_DATE2DOY DATEVEC PVL_LEAPYEAR

p = inputParser;
p.addRequired('year', @(x) (all((isnumeric(x))) & isvector(x) & ~isscalar(x)) | all(isscalar(x) & isnumeric(x) & ~isnan(x)));
p.addRequired('dayofyear',@(x) (all((isnan(x) | (isnumeric(x) & x >= 1)) & isvector(x))));
p.parse(year, dayofyear);

% make yr and doy into column vectors
year = year(:);
dayofyear = dayofyear(:);

% ensure that LY is a vector of the same size as doy
if (isscalar(year))
    yr = ones(size(dayofyear))*floor(year(:));
elseif numel(year) ~= numel(dayofyear)
    error('Error in pvl_doy2date. Input year must be scalar OR vector of same size as dayofyear.')
else
    yr = year;
end

% Figure out which years are leap years
LY = pvl_leapyear(yr);

% Determine the maximum number of days in any given year.
maxdaysinyear = 365 .* ones(size(dayofyear)) + 1 .* LY;

% Check to make sure that all days of the year are valid for the given
% year(s)
if any(dayofyear >= (maxdaysinyear+1))
    error('Error in pvl_doy2date. Day of year may not be >=366 (>=367 for leap years).')
end

% Create a matrix of size [numel(doy) x 12]. Where each row corresponds to
% the day of the year of the start of each month, accounting for leap
% years.
monthstartdays = ones(numel(dayofyear),1) * [1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335];
monthstartdays = monthstartdays + LY*[0 0 1 1 1 1 1 1 1 1 1 1];

% Multiply the doy vector with a vector of ones to create an array of size
% [numel(doy) x 12], where all columns of each row contain the DOY. Then
% compare that with the array of starting days. For each DOY value greater
% than the starting days value, it is in that month or a higher month.
mo=sum(((dayofyear * ones(1,12)) >= monthstartdays) , 2);
mo(isnan(dayofyear))=1; % If any input dayofyear is NaN, set mo to be nonzero (for indexing later)

% Subtract the doys associated with each month from the doys and add one to
% get the Day of Month (da).
MonthOffset = [1;2;0;1;1;2;2;3;4;4;5;5]; %Offset from (n-1)*30
MonthStart = ((mo - 1) .* 30) + MonthOffset(mo)+ (mo>2).*LY;
da = 1 + dayofyear-MonthStart;
da(isnan(dayofyear) | isnan(year))=NaN;
mo(isnan(dayofyear) | isnan(year))=NaN;
yr(isnan(dayofyear) | isnan(year))=NaN;