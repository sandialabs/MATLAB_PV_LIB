function TMYData = pvl_readtmy3(varargin)
% PVL_READTMY3 Read a TMY3 file in to a MATLAB struct
%
% Syntax
%   TMYData = pvl_readtmy3()
%   TMYData = pvl_readtmy3(FileName)
%
% Description
%   Read a TMY3 file and make a struct of the data. Note that values
%   contained in the struct are unchanged from the TMY3 file (i.e. units 
%   are retained). In the case of any discrepencies between this
%   documentation and the TMY3 User's Manual ([1]), the TMY3 User's Manual
%   takes precedence.
%
%   If a FileName is not provided, the user will be prompted to browse to
%   an appropriate TMY3 file.
%
%   Input
%     FileName - an optional argument which allows the user to select which
%     TMY3 format file should be read. A file path may also be necessary if
%     the desired TMY3 file is not in the MATLAB working path.
%
%   Output
%     A struct, TMYData, is provided with the following components. Note
%     that for more detailed descriptions of each component, please consult
%     the TMY3 User's Manual ([1]), especially tables 1-1 through 1-6. If
%     the output size is not specified, it is a 8760x1 vector of type double.
%
%       TMYData.SiteID - Site identifier code (USAF number), scalar double
%       TMYData.StationName - Station name, 1x1 cell string
%       TMYData.StationState - Station state 2 letter designator, 1x1 cell string
%       TMYData.SiteTimeZone - Hours from Greenwich, scalar double
%       TMYData.SiteLatitude - Latitude in decimal degrees, scalar double
%       TMYData.SiteLongitude - Longitude in decimal degrees, scalar double
%       TMYData.SiteElevation - Site elevation in meters, scalar double
%       TMYData.DateString - Date string in mm/dd/yyyy format, 8760x1 cell string
%       TMYData.TimeString - Time string in HH:MM format, local standard time, 8760x1 cell string
%       TMYData.DateNumber - Combination of date/time in MATLAB serial date (datenum) format, 8760x1 double
%       TMYData.ETR - Extraterrestrial horizontal radiation recv'd during 60 minutes prior to timestamp, Wh/m^2
%       TMYData.ETRN - Extraterrestrial normal radiation recv'd during 60 minutes prior to timestamp, Wh/m^2
%       TMYData.GHI - Direct and diffuse horizontal radiation recv'd during 60 minutes prior to timestamp, Wh/m^2
%       TMYData.GHISource - See [1], Table 1-4
%       TMYData.GHIUncertainty - Uncertainty based on random and bias error estimates - see [2]
%       TMYData.DNI - Amount of direct normal radiation (modeled) recv'd during 60 mintues prior to timestamp, Wh/m^2
%       TMYData.DNISource - See [1], Table 1-4
%       TMYData.DNIUncertainty - Uncertainty based on random and bias error estimates - see [2]
%       TMYData.DHI - Amount of diffuse horizontal radiation recv'd during 60 minutes prior to timestamp, Wh/m^2
%       TMYData.DHISource - See [1], Table 1-4
%       TMYData.DHIUncertainty - Uncertainty based on random and bias error estimates - see [2]
%       TMYData.GHillum - Avg. total horizontal illuminance recv'd during the 60 minutes prior to timestamp, lx
%       TMYData.GHillumSource - See [1], Table 1-4
%       TMYData.GHillumUncertainty - Uncertainty based on random and bias error estimates - see [2]
%       TMYData.DNillum - Avg. direct normal illuminance recv'd during the 60 minutes prior to timestamp, lx
%       TMYData.DNillumSource - See [1], Table 1-4
%       TMYData.DNillumUncertainty - Uncertainty based on random and bias error estimates - see [2]
%       TMYData.DHillum - Avg. horizontal diffuse illuminance recv'd during the 60 minutes prior to timestamp, lx
%       TMYData.DHillumSource - See [1], Table 1-4
%       TMYData.DHillumUncertainty - Uncertainty based on random and bias error estimates - see [2]
%       TMYData.Zenithlum - Avg. luminance at the sky's zenith during the 60 minutes prior to timestamp, cd/m^2
%       TMYData.ZenithlumSource - See [1], Table 1-4
%       TMYData.ZenithlumUncertainty - Uncertainty based on random and bias error estimates - see [1] section 2.10
%       TMYData.TotCld - Amount of sky dome covered by clouds or obscuring phenonema at time stamp, tenths of sky
%       TMYData.TotCldSource - See [1], Table 1-5, 8760x1 cell array of strings
%       TMYData.TotCldUnertainty - See [1], Table 1-6
%       TMYData.OpqCld - Amount of sky dome covered by clouds or obscuring phenonema that prevent observing the 
%         sky at time stamp, tenths of sky
%       TMYData.OpqCldSource - See [1], Table 1-5, 8760x1 cell array of strings
%       TMYData.OpqCldUncertainty - See [1], Table 1-6
%       TMYData.DryBulb - Dry bulb temperature at the time indicated, deg C
%       TMYData.DryBulbSource - See [1], Table 1-5, 8760x1 cell array of strings
%       TMYData.DryBulbUncertainty - See [1], Table 1-6
%       TMYData.DewPoint - Dew-point temperature at the time indicated, deg C
%       TMYData.DewPointSource - See [1], Table 1-5, 8760x1 cell array of strings
%       TMYData.DewPointUncertainty - See [1], Table 1-6
%       TMYData.RHum - Relative humidity at the time indicated, percent
%       TMYData.RHumSource - See [1], Table 1-5, 8760x1 cell array of strings
%       TMYData.RHumUncertainty - See [1], Table 1-6
%       TMYData.Pressure - Station pressure at the time indicated, 1 mbar
%       TMYData.PressureSource - See [1], Table 1-5, 8760x1 cell array of strings
%       TMYData.PressureUncertainty - See [1], Table 1-6
%       TMYData.Wdir - Wind direction at time indicated, degrees from north (360 = north; 0 = undefined,calm) 
%       TMYData.WdirSource - See [1], Table 1-5, 8760x1 cell array of strings
%       TMYData.WdirUncertainty - See [1], Table 1-6
%       TMYData.Wspd - Wind speed at the time indicated, meter/second
%       TMYData.WspdSource - See [1], Table 1-5, 8760x1 cell array of strings
%       TMYData.WspdUncertainty - See [1], Table 1-6
%       TMYData.Hvis - Distance to discernable remote objects at time indicated (7777=unlimited), meter
%       TMYData.HvisSource - See [1], Table 1-5, 8760x1 cell array of strings
%       TMYData.HvisUncertainty - See [1], Table 1-6
%       TMYData.CeilHgt - Height of cloud base above local terrain (7777=unlimited), meter
%       TMYData.CeilHgtSource - See [1], Table 1-5, 8760x1 cell array of strings
%       TMYData.CeilHgtUncertainty - See [1], Table 1-6
%       TMYData.Pwat - Total precipitable water contained in a column of unit cross section from 
%         earth to top of atmosphere, cm
%       TMYData.PwatSource - See [1], Table 1-5, 8760x1 cell array of strings
%       TMYData.PwatUncertainty - See [1], Table 1-6
%       TMYData.AOD - The broadband aerosol optical depth per unit of air mass due to extinction by
%         aerosol component of atmosphere, unitless
%       TMYData.AODSource - See [1], Table 1-5, 8760x1 cell array of strings
%       TMYData.AODUncertainty - See [1], Table 1-6
%       TMYData.Alb - The ratio of reflected solar irradiance to global horizontal irradiance, unitless
%       TMYData.AlbSource - See [1], Table 1-5, 8760x1 cell array of strings
%       TMYData.AlbUncertainty - See [1], Table 1-6
%       TMYData.Lprecipdepth - The amount of liquid precipitation observed at indicated time for the period indicated 
%         in the liquid precipitation quantity field, millimeter
%       TMYData.Lprecipquantity - The period of accumulation for the liquid precipitation depth field, hour
%       TMYData.LprecipSource - See [1], Table 1-5, 8760x1 cell array of strings
%       TMYData.LprecipUncertainty - See [1], Table 1-6
%
% Reference
% [1] Wilcox, S and Marion, W., 2008. Users Manual for TMY3 Data Sets, 
% NREL/TP-581-43156, National Renewable Energy Laboratory.  Available at
% <http://www.nrel.gov/docs/fy08osti/43156.pdf>.
%
%
% See also
%   DATEVEC  PVL_MAKELOCATIONSTRUCT  PVL_MAKETIMESTRUCT  PVL_READTMY2


if (size(varargin) == 0)
    [FileNameAndExt, FilePath]=uigetfile({ '*.csv;*.tmy3' , 'TMY3 Files (*.csv, *.tmy3)';'*.*', 'All Files (*.*)'}, 'Select a TMY3 File');
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
Header1 = textscan(FileID, '%f%q%q%f%f%f%f', 1,'Delimiter',',');
Header2 = textscan(FileID, '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s',1,'Delimiter',',');
DataLines = textscan(FileID, '%s%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%s%f%f%s%f%f%s%f%f%s%f%f%s%f%f%s%f%f%s%f%f%s%f%f%s%f%f%s%f%f%s%f%f%s%f%f%s%f%f%f%s%f',8760,'Delimiter',',');
ST = fclose(FileID);
%% Read in the data from the struct created from the textscan
TMYData.SiteID = Header1{1};
TMYData.StationName = Header1{2};
TMYData.StationState = Header1{3};
TMYData.SiteTimeZone = Header1{4};
TMYData.SiteLatitude = Header1{5};
TMYData.SiteLongitude = Header1{6};
TMYData.SiteElevation = Header1{7};

TMYData.DateString = DataLines{1};
TMYData.TimeString = DataLines{2};
TMYData.ETR = DataLines{3};
TMYData.ETRN = DataLines{4};
TMYData.GHI = DataLines{5};
TMYData.GHISource = DataLines{6};
TMYData.GHIUncertainty = DataLines{7};
TMYData.DNI = DataLines{8};
TMYData.DNISource = DataLines{9};
TMYData.DNIUncertainty = DataLines{10};
TMYData.DHI = DataLines{11};
TMYData.DHISource = DataLines{12};
TMYData.DHIUncertainty = DataLines{13};
TMYData.GHillum = DataLines{14};
TMYData.GHillumSource = DataLines{15};
TMYData.GHillumUncertainty = DataLines{16};
TMYData.DNillum = DataLines{17};
TMYData.DNillumSource = DataLines{18};
TMYData.DNillumUncertainty = DataLines{19};
TMYData.DHillum = DataLines{20};
TMYData.DHillumSource = DataLines{21};
TMYData.DHillumUncertainty = DataLines{22};
TMYData.Zenithlum = DataLines{23};
TMYData.ZenithlumSource = DataLines{24};
TMYData.ZenithlumUncertainty = DataLines{25};
TMYData.TotCld = DataLines{26};
TMYData.TotCldSource = DataLines{27};
TMYData.TotCldUnertainty = DataLines{28};
TMYData.OpqCld = DataLines{29};
TMYData.OpqCldSource = DataLines{30};
TMYData.OpqCldUncertainty = DataLines{31};
TMYData.DryBulb = DataLines{32};
TMYData.DryBulbSource = DataLines{33};
TMYData.DryBulbUncertainty = DataLines{34};
TMYData.DewPoint = DataLines{35};
TMYData.DewPointSource = DataLines{36};
TMYData.DewPointUncertainty = DataLines{37};
TMYData.RHum = DataLines{38};
TMYData.RHumSource = DataLines{39};
TMYData.RHumUncertainty = DataLines{40};
TMYData.Pressure = DataLines{41};
TMYData.PressureSource = DataLines{42};
TMYData.PressureUncertainty = DataLines{43};
TMYData.Wdir = DataLines{44};
TMYData.WdirSource = DataLines{45};
TMYData.WdirUncertainty = DataLines{46};
TMYData.Wspd = DataLines{47};
TMYData.WspdSource = DataLines{48};
TMYData.WspdUncertainty = DataLines{49};
TMYData.Hvis = DataLines{50};
TMYData.HvisSource = DataLines{51};
TMYData.HvisUncertainty = DataLines{52};
TMYData.CeilHgt = DataLines{53};
TMYData.CeilHgtSource = DataLines{54};
TMYData.CeilHgtUncertainty = DataLines{55};
TMYData.Pwat = DataLines{56};
TMYData.PwatSource = DataLines{57};
TMYData.PwatUncertainty = DataLines{58};
TMYData.AOD = DataLines{59};
TMYData.AODSource = DataLines{60};
TMYData.AODUncertainty = DataLines{61};
TMYData.Alb = DataLines{62};
TMYData.AlbSource = DataLines{63};
TMYData.AlbUncertainty = DataLines{64};
TMYData.Lprecipdepth = DataLines{65};
TMYData.Lprecipquantity = DataLines{66};
TMYData.LprecipSource = DataLines{67};
TMYData.LprecipUncertainty = DataLines{68};
%% Create a MATLAB datenum  from the string date and time
TMYData.DateNumber = datenum(strcat(TMYData.DateString,'-',TMYData.TimeString),'mm/dd/yyyy-HH:MM'); %MATLAB datenum format