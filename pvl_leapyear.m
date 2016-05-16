function LY=pvl_leapyear(Year)
% PVL_LEAPYEAR Determine if a given year is a leap year using 400 year cycles
%
% Syntax
%   LY=pvl_leapyear(Year)
%
% Description
%   The leapyear function determines if a given year, Year, in the Gregorian
%   calendar (also known as the Western calendar or Christian calendar) is a
%   leap year (a year with 29 days in the second month).
%
%   Inputs:
%       "Year" must be numeric, but may be of any size. If Year is not an
%       integer, the fractional part will be removed (floor function).
%   
%   Output:
%       "LY" is a boolean true 1 if the year is a leap year, and a
%       boolean false if the year is not a leap year. LY is of the same size as
%       Year, with LY(i) corresponding to Year(i).
%
% See also PVL_DATE2DOY PVL_DOY2DATE 
p = inputParser;
p.addRequired('Year',@(x) all(isnumeric(x)));
p.parse(Year);  

Year = floor(Year);
LY = ((mod(Year,4)==0 & mod(Year,100)~= 0) | mod(Year,400) == 0);