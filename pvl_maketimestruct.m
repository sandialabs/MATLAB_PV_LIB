function Time = pvl_maketimestruct(datenum, UTC)
% PVL_MAKETIMESTRUCT Generate a time structure from MATLAB datenum and UTC offset code
%
% Syntax
%   Time = pvl_maketimestruct(datenum, UTC)
%
% Description
%   Generate a time structure given a MATLAB datenum and UTC offset code.
%   The time structure is used in some functions such as PVL_EPHEMERIS. The
%   time struct pairs a set of times with a UTC Offset code, such that any
%   time can be corrected to UTC.
% 
%   The input "datenum" is a vector of MATLAB serial date numbers. 
%   The input "UTC" is either a scalar corresponding to ALL of the input
%   date numbers, OR it must be a vector of the same size as "datenum" with
%   each value of "UTC" corresponding to a value of "datenum". Values of
%   UTC must be greater than or equal to -12 (Yankee Time), and must be 
%   less than or equal to 13.75 (Chatham Island Daylight Time). Through the
%   use of a vector of UTC codes one can account for daylight saving time
%   by adjusting the UTC offsets associated with values in daylight saving.
%   Note that the UTC offset convention is for positive UTC offsets for
%   locations east of the prime meridian (e.g. Eastern Standard Time is  
%   UTC -5)
%
%   The output "Time" is a structure with the following elements, note that
%   all elements of the structure are column vectors, including UTCOffset.
%
%     Time.year = The year in the gregorian calendar
%     Time.month = the month of the year (January = 1 to December = 12)
%     Time.day = the day of the month
%     Time.hour = the hour of the day
%     Time.minute = the minute of the hour
%     Time.second = the second of the minute
%     Time.UTCOffset = the UTC offset code. Using the convention that a
%       positive UTC offset is for time zones east of the prime meridian
%       (e.g. EST = -5)
%
% See also DATENUM PVL_EPHEMERIS PVL_EXCELTIME2MATLAB DATEVEC

p = inputParser;
p.addRequired('datenum',@(x) all(isnumeric(x) & isvector(x)));
p.addRequired('UTC', @(x) all(isnumeric(x) & (x<=13.75 & x>-12) & isvector(x)));
p.parse(datenum,UTC);

datenum = datenum(:);

if isscalar(UTC)
    Time.UTCOffset = UTC*ones(numel(datenum),1);
elseif (numel(UTC) ~= numel(datenum))
    error('Error in pvl_maketimestruct. UTC is neither a scalar nor is it the same size as datenum.')
else
    Time.UTCOffset = UTC(:);
end

[Time.year, Time.month, Time.day, Time.hour, Time.minute, Time.second] = datevec(datenum);

end