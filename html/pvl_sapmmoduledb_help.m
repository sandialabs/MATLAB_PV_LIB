%% pvl_sapmmoduledb
% Retrieves a set of SAPM coefficients from an Excel spreadsheet.
%
%% Syntax
%%
% |SAPMparam = pvl_sapmmoduledb(Entry)|
%
% |SAPMparam = pvl_sapmmoduledb(Entry, DBfile)|
%
%% Description
% Retrieves a set of Sandia Array Performance Model (SAPM) coefficients
% from an Excel spreadsheet. The workbook must be in the same format 
% as the "SandiaModuleDatabase_20111114.xlsx" example file included with 
% the PVLib toolbox. A description of each of the coefficients associated
% with SAPM is given in the Output section of this help document, along
% with the column numbers of the data within the Excel sheet.
%
%% Inputs
%%
% * *|Entry|* - The entry number of the desired module. Entry is the numeric
%     order of the PV modules in the spreadsheet, not including header rows.
%     For example, if there is a single row of text-based header information
%     in the spreadsheet, and the user desires the module in row number 35,
%     then Entry = 34 = (35 - 1) where the "-1" accounts for the single
%     text-based header row.
% * *|DBfile|* - an optional argument which allows the user to select which
%     "database" file should be read. A file path may also be necessary if
%     the desired file is not in the MATLAB working path. If the input
%     DBfile is not included, the user will be prompted to select a file
%     using a browsing window. The file should be an Excel workbook with
%     the module data in a sheet called "Data". Row 1 should be text
%     headers. Rows 2-n contain the data specified below in columns 1
%     through 42 (A through AQ).
%
%% Outputs
%%
% * *|SAPMparam|* - A struct of SAPM module performance parameters [1]. Note that the reference
%     condition specified for most modules is Ee = 1, cell temperature = 25C, 
%     AMa = 1.5.
%%
% * *|SAPMparam.name|* = Module model name and manufacturer. Given in column 1
% * *|SAPMparam.vintage|* = Module vintage, approximate year of manufacture or
%     test. '(E)' denotes modules which have estimated parameters from
%     testing of a similar module. Given in Column 2.
% * *|SAPMparam.material|* = PV semiconductor material text found in column
% 4.
% * *|SAPMparam.area|* = Area of module in m^2 measured to outside of frame. 
%     Numeric found in column 3.
% * *|SAPMparam.AlphaIsc|* = Temperature coefficient of Isc in units of (1/C) 
%     at reference temperature (usually 25C). Numeric found in column 11.
% * *|SAPMparam.AlphaImp|* = Temperature coefficient of Imp in units of (1/C) 
%     at reference temperature (usually 25C). Numeric found in column 12.
% * *|SAPMparam.Isc0|* = Short circuit current at reference condition in units
%     of amperes. Numeric found in column 7.
% * *|SAPMparam.Imp0|* = Maximum power current at reference condition in units
%     of amperes. Numeric found in column 9.
% * *|SAPMparam.Voc0|* = Open circuit voltage at reference condition in units
%     of volts. Numeric found in column 8.
% * *|SAPMparam.Vmp0|* = Maximum power voltage at reference condition in units
%     of volts. Numeric found in column 10.
% * *|SAPMparam.BetaVoc|* = Open circuit voltage temperature coefficient at
%     reference condition (V/C). Numeric found in column 15.
% * *|SAPMparam.BetaVmp|* = Maximum power voltage temperature coefficient at
%     reference condition (V/C). Numeric found in column 17.
% * *|SAPMparam.mBetaVoc|* = Coefficient providing the irradiance dependence for
%     the BetaVoc temperature coefficient at reference irradiance (V/C).  
%     Numeric found in column 16.
% * *|SAPMparam.mBetaVmp|* = Coefficient providing the irradiance dependence for
%     the BetaVmp temperature coefficient at reference irradiance (V/C).  
%     Numeric found in column 17.
% * *|SAPMparam.Ns|* = Number of cells in series in a module's cell string(s),
%     numeric found in column 5
% * *|SAPMparam.Np|* = Number of parallel cell strings in the module, numeric
%     found in column 6
% * *|SAPMparam.delT|* = Temperature difference, in C, between module's
%     backside and the cells at reference irradiance. Typically more 
%     associated with the array mounting configuration than the module's 
%     construction properties. Numeric found in column 33.
% * *|SAPMparam.fd|* = Fraction of diffuse irradiance used by the module
%     (unitless). Typically assumed to be 1 for flat-plate modules. For
%     point-focus concentrator modules, a value of 0 is assumed, and for
%     low-concentration modules a value between 0 and 1 can be determined.
%     Numeric found in column 34.
% * *|SAPMparam.n|* = Empirically determined "diode factor" (dimensionless).
%     Numeric found in column 19.
% * *|SAPMparam.Ix0|* = Current, in amperes, of the point at 0.5*Voc at
%     reference conditions. Numeric found in column 39.
% * *|SAPMparam.Ixx0|* = Current, in amperes, of the point at 0.5*(Voc+Vmp) at
%     reference conditions. Numeric found in column 40.
% * *|SAPMparam.a_wind|* = parameter for establishing the upper limit for module 
%     temperature at low wind speeds and high solar irradiance (unitless).
%     This coefficient is strongly influenced by the array mounting 
%     configuration, and should be adjusted accordingly. Numeric found in 
%     column 35.
% * *|SAPMparam.b_wind|* = parameter for establishing the rate at which the  
%     module temperature drops as wind speed increases (unitless).
%     This coefficient is strongly influenced by the array mounting 
%     configuration, and should be adjusted accordingly. Numeric found in 
%     column 36.
%%
% * *|SAPMparam.c|* = a 1x8 vector of C values described in [1].
%%
% * *|SAPMparam.c(1)|* = C0, in column 13, value relating Imp to Ee
% (unitless).
% * *|SAPMparam.c(2)|* = C1, in column 14, value relating Imp to Ee
% (unitless).
% * *|SAPMparam.c(3)|* = C2, in column 20, value relating Vmp to Ee (unitless).
% * *|SAPMparam.c(4)|* = C3, in column 21, value relating Vmp to Ee (1/V).
% * *|SAPMparam.c(5)|* = C4, in column 37, value relating Ix to Ee
% (unitless).
% * *|SAPMparam.c(6)|* = C5, in column 38, value relating Ix to Ee
% (unitless).
% * *|SAPMparam.c(7)|* = C6, in column 41, value relating Ixx to Ee
% (unitless).
% * *|SAPMparam.c(8)|* = C7, in column 42, value relating Ixx to Ee
% (unitless).
%%
% * *|SAPMparam.a|* = a 1x5 vector of polynomial coefficients relating module 
%     current response to absolute (pressure-corrected) airmass.
% * *|SAPMparam.a(5)|* = A0, in column 22, constant coefficient of f1(AMa).
% * *|SAPMparam.a(4)|* = A1, in column 23, 1st order coefficient of
% f1(AMa).
% * *|SAPMparam.a(3)|* = A2, in column 24, 2nd order coefficient of
% f1(AMa).
% * *|SAPMparam.a(2)|* = A3, in column 25, 3rd order coefficient of
% f1(AMa).
% * *|SAPMparam.a(1)|* = A4, in column 26, 4th order coefficient of
% f1(AMa).
%%
% * *|SAPMparam.b|* = a 1x6 vector of polynomial coefficients relating module 
%     current response to incident angle (AOI).
%%
% * *|SAPMparam.b(6)|* = B0, in column 27, constant coefficient of f2(AOI).
% * *|SAPMparam.b(5)|* = B1, in column 28, 1st order coefficient of
% f2(AOI).
% * *|SAPMparam.b(4)|* = B2, in column 29, 2nd order coefficient of
% f2(AOI).
% * *|SAPMparam.b(3)|* = B3, in column 30, 3rd order coefficient of
% f2(AOI).
% * *|SAPMparam.b(2)|* = B4, in column 31, 4th order coefficient of
% f2(AOI).
% * *|SAPMparam.b(1)|* = B5, in column 32, 5th order coefficient of
% f2(AOI).
%
%% Example
%
DBfile = 'SandiaModuleDatabase_20120925.xlsx';
SAPMparam = pvl_sapmmoduledb(123, DBfile)
%% References
%%
% [1] King, D. et al, 2004. Sandia Photovoltaic Array Performance Model, SAND2004-3535,
% Sandia National Laboratories, Albuquerque, NM.
% <http://energy.sandia.gov/wp/wp-content/gallery/uploads/SAND-2004_PV-Performance-Array-Model.pdf Web Link>
%% See Also 
% <pvl_sapm_help.html |pvl_sapm|>,
% <pvl_sapmcelltemp_help.html |pvl_sapmcelltemp|>
%%
% Copyright 2014 Sandia National Laboratories
