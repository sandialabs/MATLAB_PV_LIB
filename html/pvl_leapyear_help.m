%% pvl_leapyear 
% Determine if a given year is a leap year using 400 year cycles.
%
%% Syntax
% |LY = pvl_leapyear(Year)| 
%
%% Description
% Determines if a given year, |Year|, in the Gregorian
% calendar is a leap year (a year with 29 days in the second month).
%
%% Inputs
%%
% * *|Year|* - must be numeric, but may be of any size. If |Year| is not an
% integer, the fractional part will be removed (floor function).
%   
%% Outputs
%%
% * *|LY|* - is a boolean true 1 if the year is a leap year, and a
% boolean false if the year is not a leap year. |LY| is of the same size as
% |Year|, with |LY(i)| corresponding to |Year(i)|.
%% Example 1
% Check whether 2012 is a leap year (is a leap year because 2012 is divisible by 4) 
pvl_leapyear(2012)
%% Example 2
% Check whether 2011 is a leap year (is not a leap year because 2011 is not divisible
% by 4)
pvl_leapyear(2011)
%% Example 3
% Check whether 1900 is a leap year (is not a leap year because 1900 is divisible by 100 but is not divisible
% by 400)
pvl_leapyear(1900)
%% See Also
% <pvl_date2doy_help.html |pvl_date2doy|>,     
% <pvl_doy2date_help.html |pvl_doy2date|>    
%%
% Copyright 2014 Sandia National Laboratories