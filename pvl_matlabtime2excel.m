function ExcTime = pvl_matlabtime2excel(MatTime)
% PVL_MATLABTIME2EXCEL Convert a MATLAB serial datenum to a time recognizable by Microsoft Excel
% 
% Syntax
%    ExcTime = pvl_matlabtime2excel(MatTime)
%
% Description
%    Converts a MATLAB serial datenum to a serial time
%    number recognizable by Microsoft Excel. Specifically, this
%    converts MATLAB's date number system, where day = 1 corresponds to 
%    Jan-1-0000 00:00:00, to the Microsoft Excel date number system where 
%    day = 1 corresponds to Jan-1-1900 00:00:00. 
%
%    The input MatTime is an array of dates in MATLAB's serial datenum format.
%
%    Output ExcTime is an array of dates in MS Excel's serial 1900 date
%    system. ExcTime is of the same size as MatTime.
%
% See also PVL_EXCELTIME2MATLAB PVL_RMBTIME2MATLAB DATENUM DATEVEC 


p = inputParser;
p.addRequired('MatTime',@(x) all(isnumeric(x)));
p.parse(MatTime);


ExcTime = MatTime - 693960;
