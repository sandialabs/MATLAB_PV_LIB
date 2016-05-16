function obst = pvl_readISH(ISHname)
% 
% pvl_readISH : read text file of Integrated Surface History (ISH) data
%
% Syntax:
%   data = readISH(ISHname)
%
% Input:
%   ISHname: string specifying file containing ISH data
%
% Output:
%   data: table containing ISH data, one row per observation,
%         one column for each mandatory data element. The
%         last column in the table contains a string with the
%         additional data and remarks.
%
%   Data columns are as follows (see [1] for code meanings):
%    Label         Type      Description (units)
%     len          int       length of variable portion of the record
%     catno        string    AF Catalog Station Number
%     wban         string    WBAN number
%     date         string    Date of observation, in YYYYMMDD format
%     time         string    Time of observation, in HHMM format
%     source       char      Data source code
%     lat          double    Observation latitude (degrees North)
%     long         double    Observation longitude (degrees East)
%     type         string    Report type code
%     elev         double    Observation elevation (meters)
%     station      string    Station call letters
%     qc           string    Name of quality control process applied
%     winddir      double    Wind direction (degrees clockwise from N)
%     winddirq     char      Wind direction quality code
%     windtype     char      Wind observation type code
%     windspeed    double    Wind speed (meters/second)
%     windspeedq   char      Wind speed quality code
%     ceil         double    Cloud ceiling height (meters AGL)
%     ceilq        char      Cloud ceiling quality code
%     ceilmethod   char      Cloud ceiling measurement method code
%     CAVOK        char      Ceiling and Visibility OK code
%     visdist      double    Visibility distance (meters)
%     visdistq     char      Visibility distance quality code
%     visvar       char      Visibility variability code
%     visvarq      char      Visibility variability quality code
%     temp         double    Air temperature (degrees C)
%     tempq        char      Air temperature quality code
%     dewpt        double    Dew point temperature (degrees C)
%     dewptq       char      Dew point temperature quality code
%     pressure     double    Barometric pressure (hectopascals)
%     pressureq    char      Barometric pressure quality
%     remainder    string    Optional part of record (ADD + REM)
%
% Reference:
%    [1] National Climatic Data Center, 14th Weather Squadron, Fleet
%    Numerical Meteorology and Oceanography Detachment, "Federal Climate 
%    Complex Data Documentation for Integrated Surface Data," 21 May 2014,
%    http://www1.ncdc.noaa.gov/pub/data/noaa/ish-format-document.pdf, accessed
%    5 Jun 2014
%
%% Note:
%    Requires Matlab 2013b or later
%
%% Notes for possible enhancements:
%  - Revise some column labels 
%  - Convert date and time fields into single MATLAB date string
%  - Parse remainder column (could be a separate function)
%  - Decode source, type, and quality codes
%  - Replace missing values with NaN
%  - Provide ability to filter based on quality codes

%% check Matlab version

vers = version;

pos = regexpi(vers,'(R');
yr = str2num(vers(pos+2:pos+5));
mod = vers(pos+6);
if yr<2013 || ((yr==2013) && mod=='a')
    display(['Matlab ' vers ' detected'])
    display('Function pvl_readISH requires R2013b or later')
    return
end

%% read ISH data file into string
fid = fopen(ISHname, 'r');
x = fread(fid,inf,'uint8=>char')';
fclose(fid);

% parse x into observations in a data structure 
% find boundaries of each observation (newline character)
idxs = [0 strfind(x, char(10))];

% number of records is number of newlines found (subtract index 0)
nobs = length(idxs)-1;

% create data structure 
dat(nobs,1)=struct('len',0, ... % length of variable portion of line
    'catno','XXXXXX',...        % AF Catalog Station Number
    'wban','XXXXX',...          % WBAN number
    'date','YYYYMMDD',...       % Date
    'time','HHMM',...           % Time
    'source','X',...            % Data source
    'lat',0.0,...               % Latitude (+ = North)
    'long',0.0,...              % Longitude (+ = East)
    'type','XXXXX',...          % Report type code
    'elev',0.0,...              % Elevation (m)
    'station','XXXXX',...       % Station call letter ID
    'qc','XXXX',...             % Quality control process name
    'winddir',0.0,...           % Wind direction (degrees clockwise from N)
    'winddirq','X',...          % Wind direction quality
    'windtype','X',...          % Wind observation type code
    'windspeed',0.0,...         % Wind speed (m/s)
    'windspeedq','X',...        % Wind speed quality
    'ceil',0.0,...              % Ceiling height (m, AGL)
    'ceilq','X',...             % Ceiling quality
    'ceilmethod','X',...        % Ceiling measurement method code
    'CAVOK','X',...             % Ceiling and Visibility OK code
    'visdist',0.0,...           % Visibility distance (m)
    'visdistq','X',...          % Visibility distance quality
    'visvar','X',...            % Visibility variability
    'visvarq','X',...           % Visibility variability quality
    'temp',0.0,...              % Air temperature (degrees C)
    'tempq','X',...             % Air temperature quality
    'dewpt',0.0,...             % Dew point temperature (degrees C)
    'dewptq','X',...            % Dew point temperature quality
    'pressure',0.0,...          % Barometric pressure (hPa)
    'pressureq','X',...         % Barometric pressure quality
    'remainder',[]);            % Remainder of record (ADD + REM)

% populate structure by reading data records observation by observation 
for ii = 1:nobs,
  s=x(idxs(ii)+1:idxs(ii+1)-1); % copy next observation 

  % parse the observation
  dat(ii).len = str2num(s(1:4));
  dat(ii).catno = s(5:10);
  dat(ii).wban = s(11:15);
  dat(ii).date = s(16:23);
  dat(ii).time = s(24:27);
  dat(ii).source = s(28);
  dat(ii).lat = str2num(s(29:34))/1000;
  dat(ii).long = str2num(s(35:41))/1000;
  dat(ii).type = s(42:46);
  dat(ii).elev = str2num(s(47:51));
  dat(ii).station = s(52:56);
  dat(ii).qc = s(57:60);
  dat(ii).winddir = str2num(s(61:63));
  dat(ii).winddirq = s(64);
  dat(ii).windtype = s(65);
  dat(ii).windspeed = str2num(s(66:69))/10;
  dat(ii).windspeedq = s(70);
  dat(ii).ceil = str2num(s(71:75));
  dat(ii).ceilq = s(76);
  dat(ii).ceilmethod = s(77);
  dat(ii).CAVOK = s(78);
  dat(ii).visdist = str2num(s(79:84));
  dat(ii).visdistq = s(85);
  dat(ii).visvar = s(86);
  dat(ii).visvarq = s(87);
  dat(ii).temp = str2num(s(88:92))/10;
  dat(ii).tempq = s(93);
  dat(ii).dewpt = str2num(s(94:98))/10;
  dat(ii).dewptq = s(99);
  dat(ii).pressure = str2num(s(100:104))/10;
  dat(ii).pressureq = s(105);
  dat(ii).remainder = s(106:end);
end  % for

% create table from structure 
obst=struct2table(dat);

% add units to variable columns 
obst.Properties.VariableUnits{'len'} = 'chars';
obst.Properties.VariableUnits{'lat'} = 'deg N';
obst.Properties.VariableUnits{'long'} = 'deg E';
obst.Properties.VariableUnits{'elev'} = 'm';
obst.Properties.VariableUnits{'winddir'} = 'deg';
obst.Properties.VariableUnits{'windspeed'} = 'm/s';
obst.Properties.VariableUnits{'ceil'} = 'm';
obst.Properties.VariableUnits{'visdist'} = 'm';
obst.Properties.VariableUnits{'temp'} = 'deg C';
obst.Properties.VariableUnits{'dewpt'} = 'deg C';
obst.Properties.VariableUnits{'pressure'} = 'hPa';
end % function


