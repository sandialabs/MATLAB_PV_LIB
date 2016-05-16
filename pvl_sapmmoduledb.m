function SAPMparam = pvl_sapmmoduledb(Entry,varargin)
%PVL_SAPMMODULEDB Retrieves Sandia PV Array Performance Model coefficients
%
% Syntax
%   SAPMparam = pvl_sapmmoduledb(Entry)
%   SAPMparam = pvl_sapmmoduledb(Entry, DBfile)
%
% Description
%   Retrieve a set of Sandia PV Array Performance Model (SAPM) coefficients
%   from an Excel spreadsheet. The workbook must be in the same format 
%   as the <SandiaModuleDatabase_20111114.xls> example file included with 
%   the toolbox. A description of each of the module parameters associated
%   with SAPM is given in the Output section of this help document, also,
%   the column numbers of the data within the Excel sheet is listed.
%
% Input
%   Entry - The entry number of the desired module. Entry is the numeric
%     order of the PV modules in the spreadsheet, not including header rows.
%     For example, if there is a single row of text-based header information
%     in the spreadsheet, and the user desires the module in row number 35,
%     then Entry = 34 = (35 - 1) where the "-1" accounts for the single
%     text-based header row.
%   DBfile - an optional argument which allows the user to select which
%     "database" file should be read. A file path may also be necessary if
%     the desired file is not in the MATLAB working path. If the input
%     DBfile is not included, the user will be prompted to select a file
%     using a browsing window. The file should be an Excel workbook with
%     the module data in a sheet called "Data". Row 1 should be text
%     headers. Rows 2-n contain the data specified below in columns 1
%     through 42 (A through AQ).
%
% Output
%   SAPMparam - A struct of SAPM module performance parameters with
%     the following components. Note that for more detailed descriptions of
%     each component, please consult reference [1]. Note that the reference
%     condition specified for most modules is Ee = 1, cell temperature = 25C, 
%     AMa = 1.5
%
%   SAPMparam.name = Module model name and manufacturer given in column 1
%   SAPMparam.vintage = Module vintage, approximate year of manufacture or
%     test. '(E)' denotes modules which have estimated parameters from
%     testing of a similar module. Given in Column 2
%   SAPMparam.material = PV semiconductor material text found in coulmn 4
%   SAPMparam.area = Area of module in m^2 measured to outside of frame,
%     numeric found in column 3.
%   SAPMparam.AlphaIsc = Temperature coefficient of Isc in units of (1/C) 
%     at reference temperature (usually 25C). Numeric found in column 11
%   SAPMparam.AlphaImp = Temperature coefficient of Imp in units of (1/C) 
%     at reference temperature (usually 25C). Numeric found in column 12
%   SAPMparam.Isc0 = Short circuit current at reference condition in units
%     of amperes. Numeric found in column 7.
%   SAPMparam.Imp0 = Maximum power current at reference condition in units
%     of amperes. Numeric found in column 9.
%   SAPMparam.Voc0 = Open circuit voltage at reference condition in units
%     of volts. Numeric found in column 8.
%   SAPMparam.Vmp0 = Maximum power voltage at reference condition in units
%     of volts. Numeric found in column 10.
%   SAPMparam.BetaVoc = Open circuit voltage temperature coefficient at
%     reference condition (V/C). Numeric found in column 15.
%   SAPMparam.BetaVmp = Maximum power voltage temperature coefficient at
%     reference condition (V/C). Numeric found in column 17.
%   SAPMparam.mBetaVoc = Coefficient providing the irradiance dependence for
%     the BetaVoc temperature coefficient at reference irradiance (V/C).  
%     Numeric found in column 16.
%   SAPMparam.mBetaVmp = Coefficient providing the irradiance dependence for
%     the BetaVmp temperature coefficient at reference irradiance (V/C).  
%     Numeric found in column 17.
%   SAPMparam.Ns = Number of cells in series in a module's cell string(s),
%     numeric found in column 5
%   SAPMparam.Np = Number of parallel cell strings in the module, numeric
%     found in column 6
%   SAPMparam.delT = Temperature difference, in C, between module's
%     backside and the cells at reference irradiance. Typically more 
%     associated with the array mounting configuration than the module's 
%     construction properties. Numeric found in column 33.
%   SAPMparam.fd = Fraction of diffuse irradiance used by the module
%     (unitless). Typically assumed to be 1 for flat-plate modules. For
%     point-focus concentrator modules, a value of 0 is assumed, and for
%     low-concentration modules a value between 0 and 1 can be determined.
%     Numeric found in column 34.
%   SAPMparam.n = Empirically determined "diode factor" (dimensionless),
%     numeric found in column 19
%   SAPMparam.Ix0 = Current, in amperes, of the point at 0.5*Voc at
%     reference conditions. Numeric found in column 39.
%   SAPMparam.Ixx0 = Current, in amperes, of the point at 0.5*(Voc+Vmp) at
%     reference conditions. Numeric found in column 40.
%   SAPMparam.a_wind = parameter for establishing the upper limit for module 
%     temperature at low wind speeds and high solar irradiance (unitless).
%     This coefficient is strongly influenced by the array mounting 
%     configuration, and should be adjusted accordingly. Numeric found in 
%     column 35.
%   SAPMparam.b_wind = parameter for establishing the rate at which the  
%     module temperature drops as wind speed increases (unitless).
%     This coefficient is strongly influenced by the array mounting 
%     configuration, and should be adjusted accordingly. Numeric found in 
%     column 36.
%   SAPMparam.c = a 1x8 vector of C values described in [1]
%     SAPMparam.c(1) = C0, in column 13, value relating Imp to Ee (unitless)
%     SAPMparam.c(2) = C1, in column 14, value relating Imp to Ee (unitless)
%     SAPMparam.c(3) = C2, in column 20, value relating Vmp to Ee (1/V)
%     SAPMparam.c(4) = C3, in column 21, value relating Vmp to Ee (1/V)
%     SAPMparam.c(5) = C4, in column 37, value relating Ix to Ee (unitless)
%     SAPMparam.c(6) = C5, in column 38, value relating Ix to Ee (unitless)
%     SAPMparam.c(7) = C6, in column 41, value relating Ixx to Ee (unitless)
%     SAPMparam.c(8) = C7, in column 42, value relating Ixx to Ee (unitless)
%   SAPMparam.a = a 1x5 vector of polynomial coefficients relating module 
%     current response to absolute (pressure-corrected) airmass.
%     SAPMparam.a(5) = A0, in column 22, constant coefficient of f1(AMa)
%     SAPMparam.a(4) = A1, in column 23, 1st order coefficient of f1(AMa)
%     SAPMparam.a(3) = A2, in column 24, 2nd order coefficient of f1(AMa)
%     SAPMparam.a(2) = A3, in column 25, 3rd order coefficient of f1(AMa)
%     SAPMparam.a(1) = A4, in column 26, 4th order coefficient of f1(AMa)
%   SAPMparam.b = a 1x6 vector of polynomial coefficients relating module 
%     current response to incident angle (AOI)
%     SAPMparam.b(6) = B0, in column 27, constant coefficient of f2(AOI)
%     SAPMparam.b(5) = B1, in column 28, 1st order coefficient of f2(AOI)
%     SAPMparam.b(4) = B2, in column 29, 2nd order coefficient of f2(AOI)
%     SAPMparam.b(3) = B3, in column 30, 3rd order coefficient of f2(AOI)
%     SAPMparam.b(2) = B4, in column 31, 4th order coefficient of f2(AOI)
%     SAPMparam.b(1) = B5, in column 32, 5th order coefficient of f2(AOI)
%
% References
%   [1] King, D. et al, 2004, "Sandia Photovoltaic Array Performance Model", SAND Report
%   3535, Sandia National Laboratories, Albuquerque, NM
%
%
% See also PVL_SAPM PVL_SAPMCELLTEMP

if (size(varargin) == 0)
    [FileNameAndExt, FilePath]=uigetfile({ '*.xls;*.xlsx' , 'Excel Files (*.xls, *.xlsx)';'*.*', 'All Files (*.*)'}, 'Select a file containing SAPM module coefficients');
    FilePathAndNameAndExt = [FilePath filesep FileNameAndExt];
else 
    FilePathAndNameAndExt = varargin{1};
end

p = inputParser;
p.addRequired('Entry', @(x) all(isscalar(x) & isnumeric(x)));
p.addRequired('FilePathAndNameAndExt', @(x) ischar(x) );
p.parse(Entry,FilePathAndNameAndExt);

% set the column offset due to the first row being entirely text. Note that
% if the second row becomes entirely text, then the offset would become 2.
offset =1;

[~,~,raw]=xlsread(FilePathAndNameAndExt,'Data',['A' num2str(Entry+offset) ':' 'AQ' num2str(Entry+offset)]);
SAPMparam.name = raw{1};
SAPMparam.vintage = raw{2};
SAPMparam.material = raw{4}; % 
SAPMparam.area = raw{3}; % in square meters
SAPMparam.AlphaIsc = raw{11}; % Tempcoefficient of Isc (1/C) at 25C
SAPMparam.AlphaImp = raw{12}; % Tempcoefficient of Imp (1/C) at 25C
SAPMparam.Isc0 = raw{7};
SAPMparam.Imp0 = raw{9};
SAPMparam.Voc0 = raw{8};
SAPMparam.Vmp0 = raw{10};
SAPMparam.BetaVoc = raw{15};
SAPMparam.BetaVmp = raw{17};
SAPMparam.mBetaVoc = raw{16};
SAPMparam.mBetaVmp = raw{18};
SAPMparam.Ns = raw{5}; % Series cells
SAPMparam.Np = raw{6}; % Parallel cell strings
SAPMparam.delT = raw{33}; % DeltaT, module to cell temperature delta
SAPMparam.fd = raw{34}; % Module's use of diffuse light
SAPMparam.n = raw{19}; % Diode Factor
SAPMparam.Ix0 = raw{39};
SAPMparam.Ixx0 = raw{40};

SAPMparam.a_wind = raw{35};
SAPMparam.b_wind = raw{36};

SAPMparam.c(1) = raw{13}; %C0
SAPMparam.c(2) = raw{14}; %C1
SAPMparam.c(3) = raw{20}; %C2
SAPMparam.c(4) = raw{21}; %C3
SAPMparam.c(5) = raw{37}; %C4
SAPMparam.c(6) = raw{38}; %C5
SAPMparam.c(7) = raw{41}; %C6
SAPMparam.c(8) = raw{42}; %C7

SAPMparam.a(5) = raw{22}; %A0
SAPMparam.a(4) = raw{23}; %A1
SAPMparam.a(3) = raw{24}; %A2
SAPMparam.a(2) = raw{25}; %A3
SAPMparam.a(1) = raw{26}; %A4

SAPMparam.b(6) = raw{27}; %B0 
SAPMparam.b(5) = raw{28}; %B1
SAPMparam.b(4) = raw{29}; %B2
SAPMparam.b(3) = raw{30}; %B3
SAPMparam.b(2) = raw{31}; %B4
SAPMparam.b(1) = raw{32}; %B5


end
