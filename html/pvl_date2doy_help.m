%% pvl_date2doy 
% Determine day of year from year, month of year, and day of month
%
%% Syntax
% |DOY = pvl_date2doy(Year, Month, Day)| 
%
%% Description
% Calculates the day of the year given the year, month of year, and day of
% month in the Gregorian calendar.
%
%% Inputs
%%
% * *|Year|* - scalar or vector of years (e.g. 1990). |Year(i)| corresponding 
%              to |Month(i)| and |Day(i)|. 
% * *|Month|* - scalar or vector of month numbers (1-12). |Month| must be >=1 and <13
% * *|Day|* - scalar or vector of day of month numbers.   
%             |Day| must be >=1 and <32. |Day| may either be a whole or fractional day 
%             (e.g. the 12th of the month at noon could be 12.5); if fractional, 
%             pvl_date2doy will maintain the fractional day.
%
%%
% |Year|, |Month|, and |Day| must be vectors of equal length. If |Year| or |Month| are 
% non-integers, the decimal portion of the values 
% will be dropped (floor function). The calculation utilizes the 400 year
% cycle for leap years.
%
%% Outputs
%%
% * *|DOY|* - A column vector of the day of year (whole or
%   fractional) as a floating point number from 1 to 366.999...
%
%%
% Note that |Day| must be >=1 and <32. The calculation does NOT check
%   to ensure that the day of month is valid for a given month. Thus
%   |pvl_date2doy(2012, 4, 31) = pvl_date2doy(2012, 5, 1)|. 
%% Example 1
% Determine day of year for January 1, 2012
pvl_date2doy(2012, 1, 1)
%% Example 2
% Determine day of year for December 31, 2012 (leap year)
pvl_date2doy(2012, 12, 31)
%% Example 3
% Determine day of year for December 31, 2011 (not a leap year)
pvl_date2doy(2011, 12, 31)


%% See Also
% <pvl_leapyear_help.html |pvl_leapyear|>, 
% <pvl_doy2date_help.html |pvl_doy2date|>
%%
% Copyright 2014 Sandia National Laboratories

