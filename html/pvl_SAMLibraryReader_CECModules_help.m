%% pvl_SAMLibraryReader_CECModules
% Open a System Advisor Model (SAM) CEC Module library
%
%% Syntax
% * |[CECModuleDB, ModuleNames] = pvl_SAMLibraryReader_CECModules()|
% * |[CECModuleDB, ModuleNames] = pvl_SAMLibraryReader_CECModules(LibraryFile)|
%
%% Description
% |pvl_SAMLibraryReader_CECModules| reads a System Advisor Model (SAM) [1]
% library of CEC Modules.
%
%
%% Input
% * *|LibraryFile|* - An optional input string to select which SAM library to
%     read. If omitted, a user dialog box will prompt the user to browse and
%     select a module library. Note that if |LibraryFile| is input to the
%     function, the standard MATLAB precedence order applies.
%
%% Output
% * *|CECModuleDB|* - The parameters in the SAM module library. CECModuleDB is
%   a column vector of length |NumberOfLibraryEntries| where each element.
%   is a struct with the following fields. These fields may be described in
%   more detail in the SAM documentation.
% * *|CECModuleDB.name|* - Name of the module entry in the library.
% * *|CECModuleDB.t_noct|* - Nominal Operating Cell Temperature (NOCT) in degrees C of
%           the module.
% * *|CECModuleDB.a_c|* - Module area in m^2.
% * *|CECModuleDB.Ns|* - Number of series cells within the module.
% * *|CECModuleDB.i_sc_ref|* - Short circuit current under reference
% conditions.
% * *|CECModuleDB.v_oc_ref|* - Open circuit voltage under reference
% conditions.
% * *|CECModuleDB.i_mp_ref|* - Maximum power current under reference
% condtions.
% * *|CECModuleDB.v_mp_ref|* - Maximum power voltage under reference
% conditions.
% * *|CECModuleDB.alpha_sc|* - Short circuit current temperature
% coefficient in A/C.
% * *|CECModuleDB.beta_oc|* - Open circuit voltage tempereature coefficient
% in V/C.
% * *|CECModuleDB.a_ref|* - modified diode ideality factor parameter at
%          reference conditions (units of eV), a_ref can be calculated from the
%          usual diode ideality factor (n), number of cells in series (Ns),
%          and cell temperature (Tcell) per equation (2) in [2].
% * *|CECModuleDB.IL_ref|* - Light-generated current (or photocurrent) 
%          in amperes at reference conditions. This value is referred to 
%          as Iph in some literature.
% * *|CECModuleDB.I0_ref|* - diode reverse saturation current in amperes, 
%          under reference conditions.
% * *|CECModuleDB.Rs_ref|* - series resistance under reference conditions (ohms)
% * *|CECModuleDB.Rsh_ref|* - shunt resistance under reference conditions (ohms)
% * *|adjust|* - Adjustment percentage to allow measured power temperature
%         coefficient to match reported power temperature coefficient [3].
% * *|CECModuleDB.gamma_r|* - Power temperature coefficient in %/C.
% * *|CECModuleDB.source|* - String describing the type of cells used in the module.
% * *|CECModuleDB.LibraryType|* - Library type as listed in the header information of
%          the SAM library file.
% * *|CECModuleDB.LibraryName|* - Library name as listed in the header information of
%          the SAM library file.
%
% * *|ModuleNames|* - A cell array of size [NumberOfLibraryEntries, 1] with the
%       names of the modules in the SAM library. 
%
%
%% Notes
% If this function detects that the selected input library is not of
% type "CECModule", it will display a warning to the command window and
% continue to try and read the library.
%
% The PV_LIB team would also like to thank the SAM team for maintaining
% the CEC module parameter library and allowing for interaction with the
% library files.
%
%% Sources
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
%% See also
% <pvl_singlediode_help.html |pvl_singlediode|>, <pvl_calcparams_CEC_help.html |pvl_calcparams|>
%
%%
% Copyright 2015 Sandia National Laboratories