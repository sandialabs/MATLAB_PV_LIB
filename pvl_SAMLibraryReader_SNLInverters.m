function [SNLInverterDB, InverterNames] = pvl_SAMLibraryReader_SNLInverters(varargin)
% PVL_SAMLIBRARYREADER_SNLINVERTERS Open a System Advisor Model (SAM) Sandia inverter library
%
% Syntax
%   [SNLInverterDB, InverterNames] = pvl_SAMLibraryReader_SNLInverters()
%   [SNLInverterDB, InverterNames] = pvl_SAMLibraryReader_SNLInverters(LibraryFile)
%
% Description
%   pvl_SAMLibraryReader_SNLInverters reads a System Advisor Model (SAM) [1]
%   library of Sandia inverter model parameters [2]. Reads from the .samlib
%   file format used by SAM2014.1.14 and earlier.
%
% Input Parameters:
%   LibraryFile - An optional input string to select which SAM library to
%     read. If omitted, a user dialog box will prompt the user to browse and
%     select a module library. Note that if LibraryFile is input to the
%     function, the standard MATLAB precedence order applies.
%
% Output:
%   SNLInverterDB - The parameters in the SAM inverter library. SNLInverterDB is
%   a column vector of size [NumberOfLibraryEntries, 1], where each element
%   is a struct with the following fields. These fields may be described in
%   more detail in [2].
%       name - Name of the module entry in the library. Generally the name
%         may consist of the manufacturer, model (and OEM if necessary), AC
%         voltage, data source, and year of entry into the library
%       Vac - Nominal AC voltage for the output of the inverter (V)
%       Pac0 - Maximum AC power "rating" for the inverter at nominal
%         operating condition (W)
%       Pdc0 - DC power level at which the AC power rating (Pac0) is
%         achieved (W)
%       Vdc0 - DC voltage level at which the AC power rating is achieved (V)
%       Ps0 - DC power required to start the inversion process, or
%         self-consumption by the inverter (W)
%       C0 - Parameter defining the curvature (parabolic) of the
%         relationship between AC power and DC power at the reference
%         operating condition (1/W)
%       C1 - Empirical coefficient allowing Pdc0 to vary linearly with DC
%         voltage input (1/V)
%       C2 - Empirical coefficient allowing Ps0 to vary linearly with DC
%         voltage input (1/V)
%       C3 - Empirical coefficient allowing C0 to vary linearly with DC
%         voltage input (1/V)
%       Pnt - Power consumed by the inverter when no DC input power is
%         present (e.g. at night), also known as "night tare" or "tare
%         loss" (W). Note that in the SAM library (and the preceding Sandia 
%         database) some Pnt values are positive, while others are 
%         negative, thus it may be necessary to use the absolute value of 
%         Pnt to ensure a consistent sign.
%       Vdcmax - Maximum DC voltage allowed to the inverter (V)
%       Idcmax - Maximum DC current allowed to the inverter (A)
%       MPPTLow - Low voltage limit for maximum power point tracking (V)
%       MPPTHi - High voltage limit for maximum power point tracking (V)
%       LibraryType - Library type as listed in the header information of
%          the SAM library file.
%       LibraryName - Library name as listed in the header information of
%          the SAM library file.
%
%   InverterNames - A cell array of size [NumberOfLibraryEntries, 1] with the
%       names of the inverters in the SAM library. 
%
%
% Notes:
%    If this function detects that the selected input library is not of
%    type "SNLInverter", it will display a warning to the command window and
%    continue to try and read the library.
%
%    The PV_LIB team would also like to thank the SAM team for maintaining
%    the inverter parameter library and allowing for interaction with the
%    library files.
%
% Sources:
%
% [1] System Advisor Model web page. https://sam.nrel.gov.
%
% [2] King, D. et al, 2007, "Performance Model for Grid-Connected 
%     Photovoltaic Inverters", SAND2007-5036, Sandia National 
%     Laboratories, Albuquerque, NM. 
%
%
% See also
%   PVL_SAMLIBRARYREADER_CECMODULES    PVL_SNLINVERTER  

%% Parse the input data
p = inputParser;
p.addOptional('LibraryFile', 0, @(x) ischar(x));
p.parse(varargin{:})

%% Enter Library structure basics
LibraryHeaderLines = 3; % Number of header lines before beginning data
LinesPerEntry = 16; % include line delimiters such as "!"
NonexistantCharacter = char(172); % A character which is NOT found within the library

% Create the empty structure which will hold all module data
emptyStruct = struct('name',[] ,  ...
    'Vac', [], 'Pac0', [], 'Pdc0', [], 'Vdc0', [], 'Ps0', [], ...
    'C0', [], 'C1', [], 'C2', [], 'C3', [], 'Pnt', [], ...
    'Vdcmax', [], 'Idcmax', [], 'MPPTLow', [], 'MPPTHi', [], ...
    'LibraryType', [], 'LibraryName', []);

defaultchecker = {'LibraryFile'};
if any(strcmp(defaultchecker,p.UsingDefaults))
    %% Ask user to get SAM library file
    [FileName, FilePath, FilterIndex] = ...
        uigetfile('*.samlib', 'Select a Sandia Inverter Library with .samlib extension.', 'MultiSelect', 'off');
    if FilterIndex ==0
        error('No .samlib file selected, exiting Library Reader')
    end
    FilePathandName = [FilePath FileName];
else
    FilePathandName = p.Results.LibraryFile;
end
%% Open the file and read in the header and data separately
FileID = fopen(FilePathandName);
HeaderDataIn = textscan(FileID, '%s', LibraryHeaderLines, 'Delimiter', NonexistantCharacter);
RawDataIn = textscan(FileID, '%s', 'Delimiter', NonexistantCharacter);%, 'HeaderLines', LibraryHeaderLines);
fclose(FileID);
HeaderDataIn = HeaderDataIn{1};
RawDataIn = RawDataIn{1};

%% Parse out the library name, library type, and number of entries from the header information
% This goes with every entry
LibraryName = textscan(char(HeaderDataIn(1)), '%*8s %s', 'Delimiter', NonexistantCharacter);
LibraryName = LibraryName{1};

% This goes with every entry
LibraryType = textscan(char(HeaderDataIn(2)), '%*s %s', 'Delimiter', ' ', 'MultipleDelimsAsOne', 1);
LibraryType = LibraryType{1};

if ~(strcmp(LibraryType, 'SandiaInverter'))
    warning(['A System Advisor Model (SAM) library which is NOT '...
        'of type "SandiaInverter" has been selected as input to ' ...
        'pvl_SAMLibraryReader_SNLInverters. Your data may not be correct.']);
end

NumberOfEntries = textscan(char(HeaderDataIn(3)), '%*s %f', 'Delimiter', ' ', 'MultipleDelimsAsOne', 1);
NumberOfEntries = NumberOfEntries{1};

% Create vector of structures which will hold the database information
InverterDatabase(1:NumberOfEntries) = emptyStruct;
InverterDatabase = InverterDatabase(:);

%% Step through each entry and extract the information from each entry into the structures
for cntr1 = 1:NumberOfEntries
    EntryRequested = cntr1;
    Base = (EntryRequested-1) * LinesPerEntry; % The line before the name of an entry is the "base" line
    
    % Extract information from the requested entry
    
    % Requires that the name field be as follows "entry NAMEOFENTRY"
    name = textscan(char(RawDataIn(Base+1)), '%*6s %s', 'Delimiter', NonexistantCharacter);
    InverterDatabase(EntryRequested).name = (name{1});
    
    % If we ever want to try and separate out the manufacturer and model
    % from the name field...
    %separatedname = textscan(char(name{1}), '%s %s', 'Delimiter', ':');
    %InverterDatabase(EntryRequested).manufacturer = separatedname{1};
    %InverterDatabase(EntryRequested).model = separatedname{2};
    
    Vac = textscan(char(RawDataIn(Base+2)), '%*s %f', 'Delimiter', '=');
    InverterDatabase(EntryRequested).Vac = Vac{1};
    
    Pac0 = textscan(char(RawDataIn(Base+3)), '%*s %f', 'Delimiter', '=');
    InverterDatabase(EntryRequested).Pac0 = Pac0{1};
    
    Pdc0 = textscan(char(RawDataIn(Base+4)), '%*s %f', 'Delimiter', '=');
    InverterDatabase(EntryRequested).Pdc0 = Pdc0{1};
    
    Vdc0 = textscan(char(RawDataIn(Base+5)), '%*s %f', 'Delimiter', '=');
    InverterDatabase(EntryRequested).Vdc0 = Vdc0{1};
    
    Ps0 = textscan(char(RawDataIn(Base+6)), '%*s %f', 'Delimiter', '=');
    InverterDatabase(EntryRequested).Ps0 = Ps0{1};
    
    C0 = textscan(char(RawDataIn(Base+7)), '%*s %f', 'Delimiter', '=');
    InverterDatabase(EntryRequested).C0 = C0{1};
    
    C1 = textscan(char(RawDataIn(Base+8)), '%*s %f', 'Delimiter', '=');
    InverterDatabase(EntryRequested).C1 = C1{1};
    
    C2 = textscan(char(RawDataIn(Base+9)), '%*s %f', 'Delimiter', '=');
    InverterDatabase(EntryRequested).C2 = C2{1};
    
    C3 = textscan(char(RawDataIn(Base+10)), '%*s %f', 'Delimiter', '=');
    InverterDatabase(EntryRequested).C3 = C3{1};
    
    Pnt = textscan(char(RawDataIn(Base+11)), '%*s %f', 'Delimiter', '=');
    InverterDatabase(EntryRequested).Pnt = Pnt{1};
    
    Vdcmax = textscan(char(RawDataIn(Base+12)), '%*s %f', 'Delimiter', '=');
    InverterDatabase(EntryRequested).Vdcmax = Vdcmax{1};
    
    Idcmax = textscan(char(RawDataIn(Base+13)), '%*s %f', 'Delimiter', '=');
    InverterDatabase(EntryRequested).Idcmax = Idcmax{1};
    
    MPPTLow = textscan(char(RawDataIn(Base+14)), '%*s %f', 'Delimiter', '=');
    InverterDatabase(EntryRequested).MPPTLow = MPPTLow{1};
    
    MPPTHi = textscan(char(RawDataIn(Base+15)), '%*s %f', 'Delimiter', '=');
    InverterDatabase(EntryRequested).MPPTHi = MPPTHi{1};
    
    InverterDatabase(EntryRequested).LibraryName = LibraryName;
    InverterDatabase(EntryRequested).LibraryType = LibraryType;
    
end

%% Write the information from the CEC Library to the output variables.
SNLInverterDB = InverterDatabase;
InverterNames = [SNLInverterDB.name]';
