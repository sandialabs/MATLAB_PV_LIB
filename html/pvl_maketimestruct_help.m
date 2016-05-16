%% pvl_maketimestruct 
% Generate a time struct from MATLAB datenum and UTC offset code.
%
%% Syntax
% |Time = pvl_maketimestruct(datenum, UTC)|
%
%% Description
% Generates a time structure given a MATLAB datenum and UTC offset code.
% The time structure is used in some PV_LIB functions such as <pvl_ephemeris_help.html |pvl_ephemeris|>. The
% time struct pairs a set of times with a UTC Offset code, such that any
% time can be converted to UTC.
% 
%% Inputs
%%
% * *|datenum|* - a vector of MATLAB serial date numbers. 
% * *|UTC|* = either a scalar corresponding to ALL of the input
% date numbers, OR it must be a vector of the same size as |datenum| with
% each value of |UTC| corresponding to a value of |datenum|. Values of
% |UTC| must be greater than or equal to -12 (Yankee Time), and must be 
% less than or equal to 13.75 (Chatham Island Daylight Time). Through the
% use of a vector of UTC codes one can account for daylight saving time
% by adjusting the UTC offsets associated with values in daylight saving.
% Note that the UTC offset convention is positive values for
% locations east of the prime meridian (e.g. Eastern Standard Time is  
% UTC -5).
%
%% Outputs
% * *|Time|* is a structure with the following elements, note that
% all elements of the structure are column vectors, including UTCOffset.
% * *|Time.year|* = The year in the gregorian calendar.
% * *|Time.month|* = the month of the year (January = 1 to December = 12).
% * *|Time.day|* = the day of the month.
% * *|Time.hour|* = the hour of the day.
% * *|Time.minute|* = the minute of the hour.
% * *|Time.second|* = the second of the minute.
% * *|Time.UTCOffset|* = the UTC offset code. Using the convention that a
%  positive UTC offset is for time zones east of the prime meridian
%  (e.g. EST = -5).
%
%% Example 1
% Create Time structure for a single date and time in Albuquerque, NM
Datenum = datenum('24-Oct-2003 12:45:07');
T = pvl_maketimestruct(Datenum,-7)
%% Example 2
% Create Time structure for a single date and time in Greenwich, UK
Datenum = datenum('29-Feb-2012 19:15:00');
T = pvl_maketimestruct(Datenum,0)
%% See Also 
% <pvl_ephemeris_help.html |pvl_ephemeris|>,  
% <pvl_exceltime2matlab.html |pvl_exceltime2matlab|> 
%%
% Copyright 2014 Sandia National Laboratories

