%% pvl_matlabtime2excel 
% Convert a MATLAB serial datenum to a Microsoft Excel serial date.
% 
%% Syntax
% |ExcTime = pvl_matlabtime2excel(MatTime)|
%
%% Description
% Converts a MATLAB serial datenum to a serial time
% number recognizable by Microsoft Excel. Specifically, this
% converts MATLAB's date number system, where day = 1 corresponds to 
% Jan-1-0000 00:00:00, to the Microsoft Excel date number system where day
% = 1 corresponds to Jan-1-1900 00:00:00. 
%
%% Inputs
%%
% * *|MatTime|* - an array of dates in MATLAB's serial datenum format.
%
%% Outputs    
%%
% * *|ExcTime|* - an array of dates in Microsoft Excel's serial 1900 date
%    system. |ExcTime| is of the same size as |MatTime|.
%
%% Example
%
MatTime = datenum('24-Oct-2003 12:45:07');
pvl_matlabtime2excel(MatTime)
%% See Also 
% <pvl_exceltime2matlab_help.html |pvl_exceltime2matlab|>  

%%
% Copyright 2014 Sandia National Laboratories

