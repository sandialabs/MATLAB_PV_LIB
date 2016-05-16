function MatTime = pvl_exceltime2matlab(ExcTime)
% PVL_EXCELTIME2MATLAB Convert a Microsoft Excel time to a MATLAB datenum
% 
% Syntax
%    ExcTime = MatTime = pvl_exceltime2matlab(ExcTime)
%
% Description
%    Converts a Microsoft Excel serial time to a MATLAB datenum.
%    Specifically, it converts from Excel's "1900" date system (days from 
%    Jan-1-1900 00:00:00)to MATLAB's datenum date system (days from
%    Jan-1-0000 00:00:00).
%
%    The input ExcTime is an array of dates in MS Excel's serial 1900 date
%    system.
%
%    Output MatTime is an array of dates in MATLAB's datenum date system.
%    MatTime is of the same size as ExcTime.
%
% See also PVL_MATLABTIME2EXCEL PVL_RMBTIME2MATLAB DATENUM DATEVEC 
p= inputParser;
p.addRequired('ExcTime',@isnumeric)
p.parse(ExcTime);

MatTime = p.Results.ExcTime + 693960;
