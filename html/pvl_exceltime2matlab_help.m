%% pvl_exceltime2matlab 
% Convert a Microsoft Excel time to a MATLAB datenum
% 
%% Syntax
% |MatTime = pvl_exceltime2matlab(ExcTime)|
%
%% Description
% Converts a Microsoft Excel serial time to a MATLAB datenum.
% Specifically, it converts from Excel's "1900" date system (days from 
% Jan-1-1900 00:00:00)to MATLAB's datenum date system (days from
% Jan-1-0000 00:00:00).
%
%% Inputs
%%
% * *|ExcTime|* - an array of dates in Microsoft Excel's serial 1900 date
%    system.
%
%% Outputs 
%%
% * *|MatTime|* - an array of dates in MATLAB's serial datenum format.
%
%% Example
% Example uses Jan 1, 1987 13:00
ExcTime = 31778.541666667;  
MatTime = pvl_exceltime2matlab(ExcTime)
datestr(MatTime)
%% See Also 
% <pvl_matlabtime2excel_help.html |pvl_matlabtime2excel|>  
%%
% Copyright 2014 Sandia National Laboratories
