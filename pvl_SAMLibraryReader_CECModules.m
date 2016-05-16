function [CECModuleDB, ModuleNames] = pvl_SAMLibraryReader_CECModules(varargin)
% PVL_SAMLIBRARYREADER_CECMODULES Open a System Advisor Model (SAM) CEC module library
%
% Syntax
%   [CECModuleDB, ModuleNames] = pvl_SAMLibraryReader_CECModules()
%   [CECModuleDB, ModuleNames] = pvl_SAMLibraryReader_CECModules(LibraryFile)
%
% Description
%   pvl_SAMLibraryReader_CECModules reads a System Advisor Model (SAM) [1]
%   library of CEC modules. CECModuleDB is a vector of structures which
%   describe each module in the SAM library, one element per module in the
%   SAM library. ModuleNames is a cell column vector of the names of each
%   module in the CEC Library, thus ModuleNames{n} is the same as
%   CECModuleDB(n).name. 
%
% Input Parameters:
%   LibraryFile - An optional input string to select which SAM library to
%     read. If omitted, a user dialog box will prompt the user to browse and
%     select a module library. Note that if LibraryFile is input to the
%     function, the standard MATLAB precedence order applies.
%
% Output:
%   CECModuleDB - The parameters in the SAM module library. CECModuleDB is
%   a column vector of size [NumberOfLibraryEntries, 1], where each element
%   is a struct with the following fields. These fields may be described in
%   more detail in the SAM documentation.
%       name - Name of the module entry in the library
%       t_noct - Nominal Operating Cell Temperature (NOCT) in degrees C of
%           the module
%       a_c - Module area in m^2
%       Ns - Number of series cells within the module
%       i_sc_ref - Short circuit current under reference conditions
%       v_oc_ref - Open circuit voltage under reference conditions
%       i_mp_ref - Maximum power current under reference condtions
%       v_mp_ref - Maximum power voltage under reference conditions
%       alpha_sc - Short circuit current temperature coefficient in A/C
%       beta_oc - Open circuit voltage tempereature coefficient in V/C
%       a_ref - modified diode ideality factor parameter at
%          reference conditions (units of eV), a_ref can be calculated from the
%          usual diode ideality factor (n), number of cells in series (Ns),
%          and cell temperature (Tcell) per equation (2) in [2].
%       IL_ref - Light-generated current (or photocurrent) 
%          in amperes at reference conditions. This value is referred to 
%          as Iph in some literature.
%       I0_ref - diode reverse saturation current in amperes, 
%          under reference conditions.
%       Rs_ref - series resistance under reference conditions (ohms)
%       Rsh_ref - shunt resistance under reference conditions (ohms)
%       adjust - Adjustment percentage to allow measured power temperature
%         coefficient to match reported power temperature coefficient [3].
%       gamma_r - Power temperature coefficient in %/C.
%       source - String describing the type of cells used in the module.
%       LibraryType - Library type as listed in the header information of
%          the SAM library file.
%       LibraryName - Library name as listed in the header information of
%          the SAM library file.
%
%   ModuleNames - A cell array of size [NumberOfLibraryEntries, 1] with the
%       names of the modules in the SAM library. 
%
%
% Notes:
%    If this function detects that the selected input library is not of
%    type "CECModule", it will display a warning to the command window and
%    continue to try and read the library.
%
%    The PV_LIB team would also like to thank the SAM team for maintaining
%    the CEC module parameter library and allowing for interaction with the
%    library files.
%
% Sources:
%
% [1] System Advisor Model web page. https://sam.nrel.gov.
%
% [2] W. De Soto et al., "Improvement and validation of a model for
%     photovoltaic array performance", Solar Energy, vol 80, pp. 78-88,
%     2006.
%
% [3] A. Dobos, "An Improved Coefficient Calculator for the California
%     Energy Commission 6 Parameter Photovoltaic Module Model", Journal of
%     Solar Energy Engineering, vol 134, 2012.
%
% See also
%   PVL_SINGLEDIODE      PVL_CALCPARAMS_DESOTO
%   PVL_SAMLIBRARYREADER_SNLINVERTERS

%% Parse the input data
p = inputParser;
p.addOptional('LibraryFile', 0, @(x) ischar(x));
p.parse(varargin{:})

%% Enter Library structure basics
LibraryHeaderLines = 3; % Number of header lines before beginning data
LinesPerEntry = 19; % include line delimiters such as "!"
NonexistantCharacter = char(172); % A character which is NOT found within the library

% Create the empty structure which will hold all module data
emptyStruct = struct('name',[] , 't_noct', [] , 'a_c', [], 'Ns', [], 'i_sc_ref', [], ...
    'v_oc_ref', [], 'i_mp_ref', [], 'v_mp_ref', [], 'alpha_sc', [], 'beta_oc', [], ...
    'a_ref', [], 'IL_ref', [], 'I0_ref', [], 'Rs_ref', [], 'Rsh_ref', [], ...
    'adjust', [], 'gamma_r', [], 'source', [], 'LibraryType', [], 'LibraryName', []);

defaultchecker = {'LibraryFile'};
if any(strcmp(defaultchecker,p.UsingDefaults))
    %% Ask user to get SAM library file
    [FileName, FilePath, FilterIndex] = ...
        uigetfile('*.samlib', 'Select a CEC Module Library with .samlib extension.', 'MultiSelect', 'off');
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

if ~(strcmp(LibraryType, 'CECModule'))
    warning(['A System Advisor Model (SAM) library which is NOT '...
        'of type "CECModule" has been selected as input to ' ...
        'pvl_SAMLibraryReader_CECModules. Your data may not be correct.']);
end

NumberOfEntries = textscan(char(HeaderDataIn(3)), '%*s %f', 'Delimiter', ' ', 'MultipleDelimsAsOne', 1);
NumberOfEntries = NumberOfEntries{1};

% Create vector of structures which will hold the database information
CECdatabase(1:NumberOfEntries) = emptyStruct;
CECdatabase = CECdatabase(:);

%% Step through each entry and extract the information from each entry into the structures
for cntr1 = 1:NumberOfEntries
    EntryRequested = cntr1;
    Base = (EntryRequested-1) * LinesPerEntry; % The line before the name of an entry is the "base" line
    
    % Extract information from the requested entry
    
    % Requires that the name field be as follows "entry NAMEOFENTRY"
    name = textscan(char(RawDataIn(Base+1)), '%*6s %s', 'Delimiter', NonexistantCharacter);
    CECdatabase(EntryRequested).name = (name{1});
    
    t_noct = textscan(char(RawDataIn(Base+2)), '%*s %f', 'Delimiter', '=');
    CECdatabase(EntryRequested).t_noct = t_noct{1};
    
    a_c = textscan(char(RawDataIn(Base+3)), '%*s %f', 'Delimiter', '=');
    CECdatabase(EntryRequested).a_c = a_c{1};
    
    n_s = textscan(char(RawDataIn(Base+4)), '%*s %f', 'Delimiter', '=');
    CECdatabase(EntryRequested).Ns = n_s{1};
    
    i_sc_ref = textscan(char(RawDataIn(Base+5)), '%*s %f', 'Delimiter', '=');
    CECdatabase(EntryRequested).i_sc_ref = i_sc_ref{1};
    
    v_oc_ref = textscan(char(RawDataIn(Base+6)), '%*s %f', 'Delimiter', '=');
    CECdatabase(EntryRequested).v_oc_ref = v_oc_ref{1};
    
    i_mp_ref = textscan(char(RawDataIn(Base+7)), '%*s %f', 'Delimiter', '=');
    CECdatabase(EntryRequested).i_mp_ref = i_mp_ref{1};
    
    v_mp_ref = textscan(char(RawDataIn(Base+8)), '%*s %f', 'Delimiter', '=');
    CECdatabase(EntryRequested).v_mp_ref = v_mp_ref{1};
    
    alpha_sc = textscan(char(RawDataIn(Base+9)), '%*s %f', 'Delimiter', '=');
    CECdatabase(EntryRequested).alpha_sc = alpha_sc{1};
    
    beta_oc = textscan(char(RawDataIn(Base+10)), '%*s %f', 'Delimiter', '=');
    CECdatabase(EntryRequested).beta_oc = beta_oc{1};
    
    a_ref = textscan(char(RawDataIn(Base+11)), '%*s %f', 'Delimiter', '=');
    CECdatabase(EntryRequested).a_ref = a_ref{1};
    
    i_l_ref = textscan(char(RawDataIn(Base+12)), '%*s %f', 'Delimiter', '=');
    CECdatabase(EntryRequested).IL_ref = i_l_ref{1};
    
    i_o_ref = textscan(char(RawDataIn(Base+13)), '%*s %f', 'Delimiter', '=');
    CECdatabase(EntryRequested).I0_ref = i_o_ref{1};
    
    r_s = textscan(char(RawDataIn(Base+14)), '%*s %f', 'Delimiter', '=');
    CECdatabase(EntryRequested).Rs_ref = r_s{1};
    
    r_sh_ref = textscan(char(RawDataIn(Base+15)), '%*s %f', 'Delimiter', '=');
    CECdatabase(EntryRequested).Rsh_ref = r_sh_ref{1};
    
    adjust = textscan(char(RawDataIn(Base+16)), '%*s %f', 'Delimiter', '=');
    CECdatabase(EntryRequested).adjust = adjust{1};
    
    gamma_r = textscan(char(RawDataIn(Base+17)), '%*s %f', 'Delimiter', '=');
    CECdatabase(EntryRequested).gamma_r = gamma_r{1};
    
    source = textscan(char(RawDataIn(Base+18)), '%*s %s', 'Delimiter', '=');
    CECdatabase(EntryRequested).source = (source{1});
    
    CECdatabase(EntryRequested).LibraryName = LibraryName;
    CECdatabase(EntryRequested).LibraryType = LibraryType;
    
end

%% Write the information from the CEC Library to the output variables.
CECModuleDB = CECdatabase;
ModuleNames = [CECModuleDB.name]';
