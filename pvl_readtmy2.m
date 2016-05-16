function TMYData = pvl_readtmy2(varargin)
% PVL_READTMY2 Read a TMY2 file in to a MATLAB struct
%
% Syntax
%   TMYData = pvl_readtmy2()
%   TMYData = pvl_readtmy2(FileName)
%
% Description
%   Read a TMY2 file and make a struct of the data. Note that values
%   contained in the struct are unchanged from the TMY2 file (i.e. units 
%   are retained). Time/Date and Location data imported from the TMY2 file
%   have been modified to a "friendlier" form conforming to modern
%   conventions (e.g. N latitude is postive, E longitude is positive, the
%   "24th" hour of any day is technically the "0th" hour of the next day).
%   In the case of any discrepencies between this documentation and the 
%   TMY2 User's Manual ([1]), the TMY2 User's Manual takes precedence.
%
%   If a FileName is not provided, the user will be prompted to browse to
%   an appropriate TMY2 file.
%
%   Input
%     FileName - an optional argument which allows the user to select which
%     TMY2 format file should be read. A file path may also be necessary if
%     the desired TMY2 file is not in the MATLAB working path. If FileName
%     is not provided, the user will be prompted to browse to an
%     appropriate TMY2 file.
%
%   Output
%     A struct, TMYData, is provided with the following components. Note
%     that for more detailed descriptions of each component, please consult
%     the TMY2 User's Manual ([1]), especially tables 3-1 through 3-6, and 
%     Appendix B. If the output size is not specified, it is an 8760x1 
%     vector of type double (float).
%
%       TMYData.SiteID - Site identifier code (WBAN number), scalar unsigned integer
%       TMYData.StationName - Station name, 1x1 cell string
%       TMYData.StationState - Station state 2 letter designator, 1x1 cell string
%       TMYData.SiteTimeZone - Hours from Greenwich, scalar double
%       TMYData.SiteLatitude - Latitude in decimal degrees, scalar double
%       TMYData.SiteLongitude - Longitude in decimal degrees, scalar double
%       TMYData.SiteElevation - Site elevation in meters, scalar double
%       TMYData.DateString - Date string in mm/dd/yy format, 8760x1 cell string
%       TMYData.TimeString - Time string in HH:MM format, local standard time, 8760x1 cell string
%       TMYData.DateNumber - Combination of date/time in MATLAB serial date (datenum) format, 8760x1 double
%       TMYData.ETR - Extraterrestrial horizontal radiation recv'd during 60 minutes prior to timestamp, Wh/m^2
%       TMYData.ETRN - Extraterrestrial normal radiation recv'd during 60 minutes prior to timestamp, Wh/m^2
%       TMYData.GHI - Direct and diffuse horizontal radiation recv'd during 60 minutes prior to timestamp, Wh/m^2
%       TMYData.GHISource - See [1], Table 3-3
%       TMYData.GHIUncertainty - See [1], Table 3-4
%       TMYData.DNI - Amount of direct normal radiation (modeled) recv'd during 60 mintues prior to timestamp, Wh/m^2
%       TMYData.DNISource - See [1], Table 3-3
%       TMYData.DNIUncertainty - See [1], Table 3-4
%       TMYData.DHI - Amount of diffuse horizontal radiation recv'd during 60 minutes prior to timestamp, Wh/m^2
%       TMYData.DHISource - See [1], Table 3-3
%       TMYData.DHIUncertainty - See [1], Table 3-4
%       TMYData.GHillum - Avg. total horizontal illuminance recv'd during
%         the 60 minutes prior to timestamp, units of 100 lux (e.g. value
%         of 50 = 5000 lux)
%       TMYData.GHillumSource - See [1], Table 3-3
%       TMYData.GHillumUncertainty - See [1], Table 3-4
%       TMYData.DNillum - Avg. direct normal illuminance recv'd during the
%         60 minutes prior to timestamp, units of 100 lux
%       TMYData.DNillumSource - See [1], Table 3-3
%       TMYData.DNillumUncertainty - See [1], Table 3-4
%       TMYData.DHillum - Avg. horizontal diffuse illuminance recv'd during
%         the 60 minutes prior to timestamp, units of 100 lux
%       TMYData.DHillumSource - See [1], Table 3-3
%       TMYData.DHillumUncertainty - See [1], Table 3-4
%       TMYData.Zenithlum - Avg. luminance at the sky's zenith during the
%         60 minutes prior to timestamp, units of 10 Cd/m^2 (e.g. value of
%         700 = 7,000 Cd/m^2)
%       TMYData.ZenithlumSource - See [1], Table 3-3
%       TMYData.ZenithlumUncertainty - See [1], Table 3-4
%       TMYData.TotCld - Amount of sky dome covered by clouds or obscuring phenonema at time stamp, tenths of sky
%       TMYData.TotCldSource - See [1], Table 3-5, 8760x1 cell array of strings
%       TMYData.TotCldUnertainty - See [1], Table 3-6
%       TMYData.OpqCld - Amount of sky dome covered by clouds or obscuring phenonema that prevent observing the 
%         sky at time stamp, tenths of sky
%       TMYData.OpqCldSource - See [1], Table 3-5, 8760x1 cell array of strings
%       TMYData.OpqCldUncertainty - See [1], Table 3-6
%       TMYData.DryBulb - Dry bulb temperature at the time indicated, in
%         tenths of degree C (e.g. 352 = 35.2 C).
%       TMYData.DryBulbSource - See [1], Table 3-5, 8760x1 cell array of strings
%       TMYData.DryBulbUncertainty - See [1], Table 3-6
%       TMYData.DewPoint - Dew-point temperature at the time indicated, in
%         tenths of degree C (e.g. 76 = 7.6 C).
%       TMYData.DewPointSource - See [1], Table 3-5, 8760x1 cell array of strings
%       TMYData.DewPointUncertainty - See [1], Table 3-6
%       TMYData.RHum - Relative humidity at the time indicated, percent
%       TMYData.RHumSource - See [1], Table 3-5, 8760x1 cell array of strings
%       TMYData.RHumUncertainty - See [1], Table 3-6
%       TMYData.Pressure - Station pressure at the time indicated, 1 mbar
%       TMYData.PressureSource - See [1], Table 3-5, 8760x1 cell array of strings
%       TMYData.PressureUncertainty - See [1], Table 3-6
%       TMYData.Wdir - Wind direction at time indicated, degrees from east
%         of north (360 = 0 = north; 90 = East; 0 = undefined,calm) 
%       TMYData.WdirSource - See [1], Table 3-5, 8760x1 cell array of strings
%       TMYData.WdirUncertainty - See [1], Table 3-6
%       TMYData.Wspd - Wind speed at the time indicated, in tenths of
%         meters/second (e.g. 212 = 21.2 m/s)
%       TMYData.WspdSource - See [1], Table 3-5, 8760x1 cell array of strings
%       TMYData.WspdUncertainty - See [1], Table 3-6
%       TMYData.Hvis - Distance to discernable remote objects at time 
%         indicated (7777=unlimited, 9999=missing data), in tenths of
%         kilometers (e.g. 341 = 34.1 km).
%       TMYData.HvisSource - See [1], Table 3-5, 8760x1 cell array of strings
%       TMYData.HvisUncertainty - See [1], Table 3-6
%       TMYData.CeilHgt - Height of cloud base above local terrain
%         (7777=unlimited, 88888=cirroform, 99999=missing data), in meters
%       TMYData.CeilHgtSource - See [1], Table 3-5, 8760x1 cell array of strings
%       TMYData.CeilHgtUncertainty - See [1], Table 3-6
%       TMYData.Pwat - Total precipitable water contained in a column of unit cross section from 
%         Earth to top of atmosphere, in millimeters
%       TMYData.PwatSource - See [1], Table 3-5, 8760x1 cell array of strings
%       TMYData.PwatUncertainty - See [1], Table 3-6
%       TMYData.AOD - The broadband aerosol optical depth (broadband
%         turbidity) in thousandths on the day indicated (e.g. 114 = 0.114)
%       TMYData.AODSource - See [1], Table 3-5, 8760x1 cell array of strings
%       TMYData.AODUncertainty - See [1], Table 3-6
%       TMYData.SnowDepth - Snow depth in centimeters on the day indicated,
%         (999 = missing data).
%       TMYData.SnowDepthSource - See [1], Table 3-5, 8760x1 cell array of strings
%       TMYData.SnowDepthUncertainty - See [1], Table 3-6
%       TMYData.LastSnowfall - Number of days since last snowfall (maximum 
%         value of 88, where 88 = 88 or greater days; 99 = missing data)
%       TMYData.LastSnowfallSource - See [1], Table 3-5, 8760x1 cell array of strings
%       TMYData.LastSnowfallUncertainty - See [1], Table 3-6
%       TMYData.PresentWeather - See [1], Appendix B, an 8760x1 cell array
%         of strings. Each string contains 10 numeric values. The string
%         can be parsed to determine each of 10 observed weather metrics.
%
% Reference
%   [1] Marion, W and Urban, K. "Wilcox, S and Marion, W. "User's Manual
%     for TMY2s". NREL 1995.
%
% See also
%   DATEVEC  PVL_MAKELOCATIONSTRUCT  PVL_MAKETIMESTRUCT  PVL_READTMY3


if (size(varargin) == 0)
    [FileNameAndExt, FilePath]=uigetfile({ '*.tm2' , 'TMY2 Files (*.tm2)';'*.*', 'All Files (*.*)'}, 'Select a TMY2 File');
    FilePathAndNameAndExt = [FilePath filesep FileNameAndExt];
elseif size(varargin) == 1
    FilePathAndNameAndExt = varargin{1};
    [FilePath, FileName, FileExt] = fileparts(FilePathAndNameAndExt);
    FileNameAndExt = [FileName FileExt];
end

p = inputParser;
p.addRequired('FilePathAndNameAndExt', @(x) ischar(x));
p.addRequired('FilePath', @(x) ischar(x));
p.parse(FilePathAndNameAndExt, FilePath)

FileID = fopen(FilePathAndNameAndExt);
Header1 = textscan(FileID, '%5u%22s%2s%3f%1s%2f%2f%1s%3f%2f%4f', 1,'Delimiter',',');
%                              1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58  59 60 61 62 63 64 65 66 67 68 69 70  
DataLines = textscan(FileID, '%2f%2f%2f%2f%4f%4f%4f%1s%1f%4f%1s%1f%4f%1s%1f%4f%1s%1f%4f%1s%1f%4f%1s%1f%4f%1s%1f%2f%1s%1f%2f%1s%1s%4f%1s%1s%4f%1s%1f%3f%1s%1f%4f%1s%1f%3f%1s%1f%3f%1s%1f%4f%1s%1f%5f%1s%1f%10s%3f%1s%1f%3f%1s%1f%3f%1s%1f%2f%1s%1f',8760);
ST = fclose(FileID);

% For some reason, the values for Dry Bulb Temperature and Dew Point
% Temperature do not read in correctly unless the preceeding values are
% read in as strings (this is probably due to the fact that the Dry Bulb
% and the dew point are the only numerical values which can be negative.
% Therefore, it is necessary to import the preceeding values as strings,
% then convert them to floating point numbers
DataLines{1,33}=str2double(DataLines{1,33});
DataLines{1,36}=str2double(DataLines{1,36});
%% Read in the data from the struct created from the textscan
TMYData.SiteID = Header1{1};
TMYData.StationName = strtrim(Header1{2});
TMYData.StationState = Header1{3};
TMYData.SiteTimeZone = Header1{4};
% Note that the string comparisons in the next few steps simply convert a 'N'
% latitude or 'E' longitude to positive numbers, while making a 'S'
% latitude or 'W' longitude negative numbers.
TMYData.SiteLatitude = (Header1{6}+Header1{7}./60)*(2.*strcmpi(Header1{5},'N')-1);
TMYData.SiteLongitude = (Header1{9}+Header1{10}./60)*(2.*strcmpi(Header1{8},'E')-1);
TMYData.SiteElevation = Header1{11};

TMYData.DateNumber = datenum([1900+DataLines{1} DataLines{2} DataLines{3} DataLines{4} 0*DataLines{1} 0*DataLines{1}]);
TMYData.DateString = cellstr(datestr(TMYData.DateNumber, 23));
TMYData.TimeString = cellstr(datestr(TMYData.DateNumber, 15));
TMYData.ETR = DataLines{5};
TMYData.ETRN = DataLines{6};
TMYData.GHI = DataLines{7};
TMYData.GHISource = DataLines{8};
TMYData.GHIUncertainty = DataLines{9};
TMYData.DNI = DataLines{10};
TMYData.DNISource = DataLines{11};
TMYData.DNIUncertainty = DataLines{12};
TMYData.DHI = DataLines{13};
TMYData.DHISource = DataLines{14};
TMYData.DHIUncertainty = DataLines{15};
TMYData.GHillum = DataLines{16};
TMYData.GHillumSource = DataLines{17};
TMYData.GHillumUncertainty = DataLines{18};
TMYData.DNillum = DataLines{19};
TMYData.DNillumSource = DataLines{20};
TMYData.DNillumUncertainty = DataLines{21};
TMYData.DHillum = DataLines{22};
TMYData.DHillumSource = DataLines{23};
TMYData.DHillumUncertainty = DataLines{24};
TMYData.Zenithlum = DataLines{25};
TMYData.ZenithlumSource = DataLines{26};
TMYData.ZenithlumUncertainty = DataLines{27};
TMYData.TotCld = DataLines{28};
TMYData.TotCldSource = DataLines{29};
TMYData.TotCldUnertainty = DataLines{30};
TMYData.OpqCld = DataLines{31};
TMYData.OpqCldSource = DataLines{32};
TMYData.OpqCldUncertainty = DataLines{33};
TMYData.DryBulb = DataLines{34};
TMYData.DryBulbSource = DataLines{35};
TMYData.DryBulbUncertainty = DataLines{36};
TMYData.DewPoint = DataLines{37};
TMYData.DewPointSource = DataLines{38};
TMYData.DewPointUncertainty = DataLines{39};
TMYData.RHum = DataLines{40};
TMYData.RHumSource = DataLines{41};
TMYData.RHumUncertainty = DataLines{42};
TMYData.Pressure = DataLines{43};
TMYData.PressureSource = DataLines{44};
TMYData.PressureUncertainty = DataLines{45};
TMYData.Wdir = DataLines{46};
TMYData.WdirSource = DataLines{47};
TMYData.WdirUncertainty = DataLines{48};
TMYData.Wspd = DataLines{49};
TMYData.WspdSource = DataLines{50};
TMYData.WspdUncertainty = DataLines{51};
TMYData.Hvis = DataLines{52};
TMYData.HvisSource = DataLines{53};
TMYData.HvisUncertainty = DataLines{54};
TMYData.CeilHgt = DataLines{55};
TMYData.CeilHgtSource = DataLines{56};
TMYData.CeilHgtUncertainty = DataLines{57};
TMYData.PresentWeather = DataLines{58};
TMYData.Pwat = DataLines{59};
TMYData.PwatSource = DataLines{60};
TMYData.PwatUncertainty = DataLines{61};
TMYData.AOD = DataLines{62};
TMYData.AODSource = DataLines{63};
TMYData.AODUncertainty = DataLines{64};
TMYData.SnowDepth = DataLines{65};
TMYData.SnowDepthSource = DataLines{66};
TMYData.SnowDepthUncertainty = DataLines{67};
TMYData.LastSnowfall = DataLines{68};
TMYData.LastSnowfallSource = DataLines{69};
TMYData.LastSnowfallUncertainty = DataLines{70};
