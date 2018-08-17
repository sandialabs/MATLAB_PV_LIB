%% pvl_SAMLibraryReader_SNLInverters
% Open a System Advisor Model (SAM) Sandia Inverter Library
%
%% Syntax
% * |[SNLInverterDB, InverterNames] = pvl_SAMLibraryReader_SNLInverters()|
% * |[SNLInverterDB, InverterNames] = pvl_SAMLibraryReader_SNLInverters(LibraryFile)|
%
%% Description
% |pvl_SAMLibraryReader_SNLInverters| reads a System Advisor Model (SAM) [1]
% library of Sandia inverter model parameters [2]. Reads from the |.samlib|
% file format used by SAM2014.1.14 and earlier.
%
%% Input
% * *|LibraryFile|* - An optional input string to select which SAM library to
%     read. If omitted, a user dialog box will prompt the user to browse and
%     select a module library. Note that if LibraryFile is input to the
%     function, the standard MATLAB precedence order applies.
%
%% Output
% * *|SNLInverterDB|* - The parameters in the SAM inverter library. SNLInverterDB is
%   a column vector of size [NumberOfLibraryEntries, 1], where each element
%   is a struct with the following fields. These fields may be described in
%   more detail in [2].
% * *|name|* - Name of the module entry in the library. Generally the name
%         may consist of the manufacturer, model (and OEM if necessary), AC
%         voltage, data source, and year of entry into the library
% * *|Vac|* - Nominal AC voltage for the output of the inverter (V)
% * *|Pac0|* - Maximum AC power "rating" for the inverter at nominal
%         operating condition (W)
% * *|Pdc0|* - DC power level at which the AC power rating (Pac0) is
%         achieved (W)
% * *|Vdc0|* - DC voltage level at which the AC power rating is achieved (V)
% * *|Ps0|* - DC power required to start the inversion process, or
%         self-consumption by the inverter (W)
% * *|C0|* - Parameter defining the curvature (parabolic) of the
%         relationship between AC power and DC power at the reference
%         operating condition (1/W)
% * *|C1|* - Empirical coefficient allowing Pdc0 to vary linearly with DC
%         voltage input (1/V)
% * *|C2|* - Empirical coefficient allowing Ps0 to vary linearly with DC
%         voltage input (1/V)
% * *|C3|* - Empirical coefficient allowing C0 to vary linearly with DC
%         voltage input (1/V)
% * *|Pnt|* - Power consumed by the inverter when no DC input power is
%         present (e.g. at night), also known as "night tare" or "tare
%         loss" (W). Note that in the SAM library (and the preceding Sandia 
%         database) some Pnt values are positive, while others are 
%         negative, thus it may be necessary to use the absolute value of 
%         Pnt to ensure a consistent sign.
% * *|Vdcmax|* - Maximum DC voltage allowed to the inverter (V)
% * *|Idcmax|* - Maximum DC current allowed to the inverter (A)
% * *|MPPTLow|* - Low voltage limit for maximum power point tracking (V)
% * *|MPPTHi|* - High voltage limit for maximum power point tracking (V)
% * *|LibraryType|* - Library type as listed in the header information of
%          the SAM library file.
% * *|LibraryName|* - Library name as listed in the header information of
%          the SAM library file.
%
% * |InverterNames|* - A cell array of size [NumberOfLibraryEntries, 1] with the
%       names of the inverters in the SAM library. 
%
%
%% Notes
% If this function detects that the selected input library is not of
% type "SNLInverter", it will display a warning to the command window and
% continue to try and read the library.
%
% The PV_LIB team would also like to thank the SAM team for maintaining
% the inverter parameter library and allowing for interaction with the
% library files.
%
%% Sources
%
% [1] System Advisor Model web page. https://sam.nrel.gov.
%
% [2] King, D. et al, 2007, "Performance Model for Grid-Connected 
%     Photovoltaic Inverters", SAND2007-5036, Sandia National 
%     Laboratories, Albuquerque, NM. 
%
%
%% See also
% <pvl_snlinverter_help.html |pvl_snlinverter|>
%
%%
% Copyright 2015 Sandia National Laboratories