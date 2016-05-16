%% pvl_doy2date
% Determines the Year, Month of year, and Day of month given the year and day of year
%
%% Syntax
% |[yr mo da] = pvl_doy2date(year, dayofyear)| 

%% Description
% Gives the year, month of year, and day of month in the Gregorian calendar, 
% when given the year and day of year. 
% 
%% Inputs
%%
% * *|year|* - must be a numeric vector or scalar of Gregorian years. 
%       If |year| is not an integer, the fractional part will be removed (floor).
% * *|dayofyear|* - must be a numeric vector with each element >=1 and
%       <366 (< 367 for leap years).
%
%% Outputs
%%
% * *|yr|* - a column vector of same size as |dayofyear| returning the 
% year input. If the input |year| was a scalar, then |yr| is a vector 
% filled with |year|.
% * *|mo|* - a column vector of month of the year (1 = January, 
% 12 = December)
% * *|da|* - is a column vector of day of the month, including any 
% fractional days given in |dayofyear|.
%
%% Example 1
% What is the date (m/d/yr) of the 60th day of 2011?
[yr mo da] = pvl_doy2date(2011, 60)

%% Example 2
% What is the date (m/d/yr) of the 60th day of 2012 (a leap year)?
[yr mo da] = pvl_doy2date(2012, 60)
%% See Also
% <pvl_date2doy_help.html |pvl_date2doy|> 
%%
% Copyright 2014 Sandia National Laboratories
