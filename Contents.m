% Toolbox SNL_PVLib provides a set of functions useful for modeling solar
% photovoltaic systems
% Version 1.3 Dec-2015


%% Time and Location Utilities
%   pvl_date2doy           - Determine day of year from year, month of year, and day of month
%   pvl_doy2date           - Determines the Year, Month of year, and Day of month given the year and day of year
%   pvl_leapyear           - Determine if a given year is a leap year using 400 year cycles
%   pvl_exceltime2matlab   - Convert a Microsoft Excel time to a MATLAB datenum
%   pvl_matlabtime2excel   - Convert a MATLAB serial datenum to a time recognizable by Microsoft Excel
%   pvl_rmbtime2matlab     - Creates a MATLAB datenum from the time convention used in Rocky Mountain Basic
%   pvl_maketimestruct     - Generate a time structure from MATLAB datenum and UTC offset code
%   pvl_makelocationstruct - Create a struct to define a site location

%% Irradiance and Atmospheric Functions
%   pvl_readtmy3           - Read a TMY3 file in to a MATLAB struct
%   pvl_readtmy2           - Read a TMY2 file in to a MATLAB struct
%   pvl_getISDdata         - Fetch data from NOAA ISD at ftp.ncdc.noaa.gov
%   pvl_readISH            - Read data fetched from ftp.ncdc.noaa.gov into a table 
%   pvl_ephemeris          - Calculates the position of the sun given time, location, and optionally pressure and temperature
%   pvl_spa                - Calculates the position of the sun given time, location, and optionally pressure and temperature
%   pvl_extraradiation     - Determine extraterrestrial radiation from day of year
%   pvl_alt2pres           - Determine site pressure from altitude
%   pvl_pres2alt           - Determine altitude from site pressure
%   pvl_relativeairmass    - Gives the relative (not pressure-corrected) airmass
%   pvl_absoluteairmass    - Determine absolute (pressure corrected) airmass from relative airmass and pressure
%   pvl_disc               - Estimate Direct Normal Irradiance from Global Horizontal Irradiance using the DISC model
%   pvl_dirint             - Determine DNI and DHI from GHI using the DIRINT modification of the DISC model
%   pvl_erbs               - Determine DNI and DHI from GHI using the Erbs model
%   pvl_louche             - Determine DNI and DHI from GHI using the Louche model
%   pvl_orgill_hollands    - Determine DNI and DHI from GHI using the Orgill and Hollands model
%   pvl_reindl_1           - Determine DNI and DHI from GHI using the Reindl_1 model
%   pvl_reindl_2           - Determine DNI and DHI from GHI using the Reindl_2 model
%   pvl_clearsky_haurwitz  - Determine clear sky GHI using the Haurwitz model
%   pvl_clearsky_ineichen  - Determine clear sky GHI using the Ineichen model
%   pvl_calcPwat           - Calculate precipitable water (cm) from ambient air temperature (C) and relatively humidity (%) 

%% Irradiance Translation Functions
%   pvl_grounddiffuse      - Estimate diffuse irradiance from ground reflections given irradiance, albedo, and surface tilt 
%   pvl_isotropicsky       - Determine diffuse irradiance from the sky on a tilted surface using isotropic sky model
%   pvl_reindl1990         - Determine diffuse irradiance from the sky on a tilted surface using Reindl's 1990 model
%   pvl_perez              - Determine diffuse irradiance from the sky on a tilted surface using one of the Perez models
%   pvl_kingdiffuse        - Determine diffuse irradiance from the sky on a tilted surface using the King model
%   pvl_klucher1979        - Determine diffuse irradiance from the sky on a tilted surface using Klucher's 1979 model
%   pvl_haydavies1980      - Determine diffuse irradiance from the sky on a tilted surface using Hay & Davies' 1980 model
%   pvl_getaoi             - Determine angle of incidence from surface tilt/azimuth and apparent sun zenith/azimuth

%% Photovoltaic System Functions
%   pvl_sapmmoduledb       - Retrieves Sandia Array Performance Model coefficients
%   pvl_SAMLibraryReader_CECModules    - Open a System Advisor Model (SAM) CEC module library
%   pvl_SAMLibraryReader_SNLInverters  - Open a System Advisor Model (SAM) v2014.1.14 or earlier inverter library
%   pvl_sapmcelltemp       - Estimate cell temperature from irradiance, windspeed, ambient temperature, and module parameters (SAPM)
%   pvl_physicaliam        - Determine the incidence angle modifier using refractive index, glazing thickness, and extinction coefficient
%   pvl_martinruiziam      - Determine the incidence angle modifier using the Martin and Ruiz incident angle model
%   pvl_ashraeiam          - Determine the incidence angle modifier using the ASHRAE incident angle model
%   pvl_FSspeccorr         - Calculate spectral mismatch modifier based on precipitable water and absolute airmass
%   pvl_calcparams_desoto  - Calculate module performance model coefficients for the De Soto single diode model
%   pvl_calcparams_CEC     - Create module performance model coefficients for the single diode model used by the CEC
%   pvl_calcparams_PVsyst  - Create module performance model coefficient for the single diode model in PVsyst version 6
%   pvl_singlediode        - Solves the single-diode equation to obtain a photovoltaic IV curve
%   pvl_sapm               - Sandia Array Performance Model to get 5 points on IV curve
%   pvl_huld               - Calculates DC power using the Huld PV module model
%   pvl_snlinverter        - Converts DC power and voltage to AC power using Sandia's Grid-Connected PV Inverter model
%   pvl_adrinverter        - Converts DC power and voltage to AC power using Anton Driesse's Grid-Connected PV Inverter model
%   pvl_singleaxis         - Determine the rotation angle of a 1 axis tracker, and sun incident angle to tracked surface 

%%  Functions for parameter estimation for PV module models
%   pvl_PVsyst_parameter_estimation - estimates parameters for the PVsyst version 6 module performance model
%   pvl_desoto_parameter_estimation - estimates parameters for the De Soto single diode module performance model
%   pvl_rectify_IV_curve            - ensures that Isc and Voc are included in a IV curve and removes duplicate voltage and current points
%   est_single_diode_param          - estimates five parameters for an IV curve using a regression on the co-content
%   calc_theta_phi_exact            - computes the arguments for the Lambert W function in the analytic solutions to the single diode equation
%   update_Io_known_n               - iterative update to the dark current (Io) value for an IV curve to better fit Voc
%   update_Rsh_fixed_pt             - iterative update to the parallel resistance (Rsh) value for an IV curve to better fit Vmp
%   Schumaker_QSpline               - fit a non-increasing, concave downward quadratic spline to IV curve data
%

%%  Numerical utilities
%   pvl_lambertw           - Compute values for the Lambert W function W(z)
%   numdiff                - Compute numerical derivatives for unequally spaced data
