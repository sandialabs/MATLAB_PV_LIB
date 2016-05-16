%% pvl_readtmy3
% Read a Typical Meteological Year 3 (TMY3) file to a MATLAB struct.
%
%% Syntax
%%
% |TMYData = pvl_readtmy3()|
%
% |TMYData = pvl_readtmy3(FileName)|
%
%% Description
% Read a TMY3 file and make a struct of the data. Note that values
% contained in the struct are unchanged from the TMY3 file (i.e. units 
% are retained). In the case of any discrepencies between this
% documentation and the TMY3 User's Manual [1], the TMY3 User's Manual
% takes precedence.
%
%
%% Inputs
%%
% * *|FileName|* - an optional argument which allows the user to select which
%     TMY3 format file should be read. A file path may also be necessary if
%     the desired TMY3 file is not in the MATLAB working path.  If |FileName|
%     is not provided, the user will be prompted to browse to
%     an appropriate TMY3 file.
%
%% Output
%
% *  *|TMYData|* - a struct with the following components. For more detailed 
%     descriptions of each component, please consult
%     the TMY3 User's Manual [1], especially tables 1-1 through 1-6.  
%     Unless otherwise stated, each field is a column vector of type double (float).
%%
% * *|TMYData.SiteID|* - Site identifier code (USAF number), scalar double.
% * *|TMYData.StationName|* - Station name, 1x1 cell string.
% * *|TMYData.StationState|* - Station state 2 letter designator, 1x1 cell
% string.
% * *|TMYData.SiteTimeZone|* - Hours from Greenwich, scalar double.
% * *|TMYData.SiteLatitude|* - Latitude in decimal degrees, scalar double.
% * *|TMYData.SiteLongitude|* - Longitude in decimal degrees, scalar
% double.
% * *|TMYData.SiteElevation|* - Site elevation in meters, scalar double.
% * *|TMYData.DateString|* - Date string in mm/dd/yy format, 8760x1 cell
% string.
% * *|TMYData.TimeString|* - Time string in HH:MM format, local standard
% time, 8760x1 cell string.
% * *|TMYData.DateNumber|* - Combination of date/time in MATLAB serial date
% (datenum) format, 8760x1 double.
% * *|TMYData.ETR|* - Extraterrestrial horizontal radiation recv'd during
% 60 minutes prior to timestamp, Wh/m^2.
% * *|TMYData.ETRN|* - Extraterrestrial normal radiation recv'd during 60
% minutes prior to timestamp, Wh/m^2.
% * *|TMYData.GHI|* - Direct and diffuse horizontal radiation recv'd during
% 60 minutes prior to timestamp, Wh/m^2.
% * *|TMYData.GHISource|* - See [1], Table 1-4.
% * *|TMYData.GHIUncertainty|* - Uncertainty based on random and bias error
% estimates - see [2].
% * *|TMYData.DNI|* - Amount of direct normal radiation (modeled) recv'd
% during 60 mintues prior to timestamp, Wh/m^2.
% * *|TMYData.DNISource|* - See [1], Table 1-4.
% * *|TMYData.DNIUncertainty|* - Uncertainty based on random and bias error
% estimates - see [2].
% * *|TMYData.DHI|* - Amount of diffuse horizontal radiation recv'd during
% 60 minutes prior to timestamp, Wh/m^2.
% * *|TMYData.DHISource|* - See [1], Table 1-4.
% * *|TMYData.DHIUncertainty|* - Uncertainty based on random and bias error
% estimates - see [2].
% * *|TMYData.GHillum|* - Average total horizontal illuminance recv'd during the 60 minutes prior to timestamp, Wh/m^2.
% * *|TMYData.GHillumSource|* - See [1], Table 1-4.
% * *|TMYData.GHillumUncertainty|* - Uncertainty based on random and bias
% error estimates - see [2].
% * *|TMYData.DNillum|* - Average direct normal illuminance recv'd during the 60 minutes prior to timestamp, Wh/m^2.
% * *|TMYData.DNillumSource|* - See [1], Table 1-4.
% * *|TMYData.DNillumUncertainty|* - Uncertainty based on random and bias
% error estimates - see [2].
% * *|TMYData.DHillum|* - Average horizontal diffuse illuminance recv'd during the 60 minutes prior to timestamp, Wh/m^2.
% * *|TMYData.DHillumSource|* - See [1], Table 1-4.
% * *|TMYData.DHillumUncertainty|* - Uncertainty based on random and bias
% error estimates - see [2].
% * *|TMYData.Zenithlum|* - Average luminance at the sky's zenith during the 60 minutes prior to timestamp, cd/m^2
% * *|TMYData.ZenithlumSource|* - See [1], Table 1-4.
% * *|TMYData.ZenithlumUncertainty|* - Uncertainty based on random and bias
% error estimates - see [1] section 2.10.
% * *|TMYData.TotCld|* - Amount of sky dome covered by clouds or obscuring
% phenonema at time stamp, tenths of sky.
% * *|TMYData.TotCldSource|* - See [1], Table 1-5, 8760x1 cell array of
% strings.
% * *|TMYData.TotCldUnertainty|* - See [1], Table 1-6.
% * *|TMYData.OpqCld|* - Amount of sky dome covered by clouds or obscuring phenonema that prevent observing the 
% * *|  sky at time stamp, tenths of sky.
% * *|TMYData.OpqCldSource|* - See [1], Table 1-5, 8760x1 cell array of
% strings.
% * *|TMYData.OpqCldUncertainty|* - See [1], Table 1-6.
% * *|TMYData.DryBulb|* - Dry bulb temperature at the time indicated, deg
% C.
% * *|TMYData.DryBulbSource|* - See [1], Table 1-5, 8760x1 cell array of
% strings.
% * *|TMYData.DryBulbUncertainty|* - See [1], Table 1-6.
% * *|TMYData.DewPoint|* - Dew-point temperature at the time indicated, deg
% C.
% * *|TMYData.DewPointSource|* - See [1], Table 1-5, 8760x1 cell array of
% strings.
% * *|TMYData.DewPointUncertainty|* - See [1], Table 1-6.
% * *|TMYData.RHum|* - Relative humidity at the time indicated, percent.
% * *|TMYData.RHumSource|* - See [1], Table 1-5, 8760x1 cell array of
% strings.
% * *|TMYData.RHumUncertainty|* - See [1], Table 1-6.
% * *|TMYData.Pressure|* - Station pressure at the time indicated, 1 mbar.
% * *|TMYData.PressureSource|* - See [1], Table 1-5, 8760x1 cell array of
% strings.
% * *|TMYData.PressureUncertainty|* - See [1], Table 1-6.
% * *|TMYData.Wdir|* - Wind direction at time indicated, degrees from north
% (360 = north; 0 = undefined,calm).
% * *|TMYData.WdirSource|* - See [1], Table 1-5, 8760x1 cell array of
% strings.
% * *|TMYData.WdirUncertainty|* - See [1], Table 1-6.
% * *|TMYData.Wspd|* - Wind speed at the time indicated, meter/second.
% * *|TMYData.WspdSource|* - See [1], Table 1-5, 8760x1 cell array of
% strings.
% * *|TMYData.WspdUncertainty|* - See [1], Table 1-6.
% * *|TMYData.Hvis|* - Distance to discernable remote objects at time
% indicated (7777=unlimited), meters.
% * *|TMYData.HvisSource|* - See [1], Table 1-5, 8760x1 cell array of
% strings.
% * *|TMYData.HvisUncertainty|* - See [1], Table 1-6.
% * *|TMYData.CeilHgt|* - Height of cloud base above local terrain
% (7777=unlimited), meters.
% * *|TMYData.CeilHgtSource|* - See [1], Table 1-5, 8760x1 cell array of
% strings.
% * *|TMYData.CeilHgtUncertainty|* - See [1], Table 1-6.
% * *|TMYData.Pwat|* - Total precipitable water contained in a column of unit cross section from 
%    earth to top of atmosphere, cm.
% * *|TMYData.PwatSource|* - See [1], Table 1-5, 8760x1 cell array of
% strings.
% * *|TMYData.PwatUncertainty|* - See [1], Table 1-6.
% * *|TMYData.AOD|* - The broadband aerosol optical depth per unit of air mass due to extinction by
%    aerosol component of atmosphere, unitless.
% * *|TMYData.AODSource|* - See [1], Table 1-5, 8760x1 cell array of
% strings.
% * *|TMYData.AODUncertainty|* - See [1], Table 1-6.
% * *|TMYData.Alb|* - The ratio of reflected solar irradiance to global
% horizontal irradiance, unitless.
% * *|TMYData.AlbSource|* - See [1], Table 1-5, 8760x1 cell array of
% strings.
% * *|TMYData.AlbUncertainty|* - See [1], Table 1-6.
% * *|TMYData.Lprecipdepth|* - The amount of liquid precipitation observed
% at indicated time for the period indicated.
%    in the liquid precipitation quantity field, millimeters.
% * *|TMYData.Lprecipquantity|* - The period of accumulation for the liquid
% precipitation depth field, hours.
% * *|TMYData.LprecipSource|* - See [1], Table 1-5, 8760x1 cell array of
% strings.
% * *|TMYData.LprecipUncertainty|* - See [1], Table 1-6.
%
%% Example
%
pvl_readtmy3('723650TY.csv')
%% References
%%
% [1] Wilcox, S and Marion, W., 2008. Users Manual for TMY3 Data Sets, 
% NREL/TP-581-43156, National Renewable Energy Laboratory.  Available at
% <http://www.nrel.gov/docs/fy08osti/43156.pdf>.
%
%% See Also
% <pvl_makelocationstruct_help.html |pvl_makelocationstruct|> ,
% <pvl_maketimestruct_help.html.html |pvl_maketimestruct|> 
%%
% Copyright 2014 Sandia National Laboratories

