%% Functions by Category
% PV_LIB Toolbox
% Version 1.32 Feb-2016
%
%
%% Time and Location Utilities
%%
% * <pvl_date2doy_help.html |pvl_date2doy|> - Gives day of year from date
% * <pvl_doy2date_help.html |pvl_doy2date|> - Gives date from day of year
% * <pvl_leapyear_help.html |pvl_leapyear|> - Boolian inticating if year is
% a leap year
% * <pvl_matlabtime2excel_help.html |pvl_matlabtime2excel|> - Converts
% matlab datetimes to excel datetime convention
% * <pvl_exceltime2matlab_help.html |pvl_exceltime2matlab|> - Converts excel
% datetimes to Matlab datetime convention
% * <pvl_maketimestruct_help.html |pvl_maketimestruct|> - Creates a |Time| structure
% * <pvl_makelocationstruct_help.html |pvl_makelocationstruct|> - Creates a |Location| structure
%
%% Irradiance and Atmospheric Functions
%%
% * <pvl_readtmy2_help.html |pvl_readtmy2|> - Read a Typical Meteological 
% Year 2 (TMY2) file in to a MATLAB struct
% * <pvl_readtmy3_help.html |pvl_readtmy3|> - Read a Typical Meteological 
% Year 3 (TMY3) file in to a MATLAB struct
% * <pvl_getSDdata_help.html |pvl_getISDdata|> - Fetch data from NOAA ISD at ftp.ncdc.noaa.gov
% * <pvl_readISH_help.html |pvl_readISH|> - Read data fetched from ftp.ncdc.noaa.gov into a table
% * <pvl_ephemeris_help.html |pvl_ephemeris|> - Position of the sun given
% date,time, and location
% * <pvl_spa_help.html |pvl_spa|> - Position of the sun given
% date,time, and location using NREL SPA function (slower but more
% accurate)
% * <pvl_extraradiation_help.html |pvl_extraradiation|> - Extraterrestrial incident radiation
% * <pvl_pres2alt_help.html |pvl_pres2alt|> - Standard altitude based
% on air pressure
% * <pvl_alt2pres_help.html |pvl_alt2pres|> - Average atmospheric pressure
% * <pvl_relativeairmass_help.html |pvl_relativeairmass|> - Relative optical airmass
% * <pvl_absoluteairmass_help.html |pvl_absoluteairmass|> - Absolute airmass
% (assumes standard pressure and site elevation)
% * <pvl_disc_help.html |pvl_disc|> - DISC model for estimating direct normal
% irradiance from global horizontal irradiance
% * <pvl_dirint_help.html |pvl_dirint|> - DIRINT adjustment to the DISC model for estimating direct normal
% irradiance from global horizontal irradiance
% * <pvl_louche_help.html |pvl_louche|> - Estimate Direct Normal Irradiance from Global Horizontal Irradiance using the LOUCHE model
% * <pvl_erbs_help.html |pvl_erbs|> - Estimate Direct Normal Irradiance from Global Horizontal Irradiance using the Erbs model
% * <pvl_orgill_hollands_help.html |pvl_orgill_hollands|> - Estimate Direct Normal Irradiance from Global Horizontal Irradiance using the Orgill and Hollands model
% * <pvl_reindl_1_help.html |pvl_reindl_1|> - Estimate Direct Normal Irradiance from Global Horizontal Irradiance using the Reindl 1 model
% * <pvl_reindl_2_help.html |pvl_reindl_2|> - Estimate Direct Normal Irradiance from Global Horizontal Irradiance using the Reindl 2 model
% * <pvl_clearsky_haurwitz_help.html |pvl_clearsky_haurwitz|> - Clear sky global horizontal irradiance model (simple)
% * <pvl_clearsky_ineichen_help.html |pvl_clearsky_ineichen|> - Clear sky irradiance (GHI, DNI, and DHI) model
% * <pvl_calcPwat_help.html |pvl_calcPwat|> - Calculate precipitable water (cm) from ambient air temperature (C) and relatively humidity (%)

%% Irradiance Translation Functions
%%
% * <pvl_grounddiffuse_help.html |pvl_grounddiffuse|> - Ground reflected diffuse irradiance
% on a tilted plane
% * <pvl_isotropicsky_help.html |pvl_isotropicsky|> - Isotropic diffuse sky irradiance 
% on a tilted surface using 
% * <pvl_reindl1990_help.html |pvl_reindl1990|> - Reindl's 1990 model of diffuse sky irradiance 
% on a tilted surface
% * <pvl_perez_help.html |pvl_perez|> - Perez's model of diffuse sky irradiance
% on a tilted surface 
% * <pvl_kingdiffuse_help.html |pvl_kingdiffuse|> - King's model of diffuse sky irradiance 
% on a tilted surface 
% * <pvl_klucher1979_help.html |pvl_klucher1979|> - Klucher's model of diffuse sky irradiance 
% on a tilted surface
% * <pvl_haydavies1980_help.html |pvl_haydavies1980|> - Hay & Davies' model of diffuse sky irradiance 
% on a tilted surface 
% * <pvl_getaoi_help.html |pvl_getaoi|> - Determine angle of incidence between tilted array surface 
% (tilt/azimuth) and apparent sun position (zenith/azimuth) 

%% Irradiance Analysis Functions
% * <pvl_detect_clear_times_help.html |pvl_detect_clear_times|> - Identify times with GHI consistent with clear sky conditions
% * <pvl_detect_shadows_help.html |pvl_detect_shadows|> - Identify shading on a GHI instrument from nearby structures such as 
% wires and poles

%% Photovoltaic System Functions
%%
% * <pvl_sapmmoduledb_help.html |pvl_sapmmoduledb|> - Retrieves coefficients for the Sandia Array Performance Model (SAPM)
% * <pvl_SAMLibraryReader_CECModules_help.html |pvl_SAMLibraryReader_CECModules|> - Open a System Advisor Model (SAM) CEC module library
% * <pvl_SAMLibraryReader_SNLInverters_help.html |pvl_SAMLibraryReader_SNLInverters|> - Open an inverter library from System Advisor Model (SAM) v2014.1.14 or earlier 
% * <pvl_sapmcelltemp_help.html |pvl_sapmcelltemp|> - Estimate cell temperature from irradiance, windspeed, ambient temperature, and module parameters
% * <pvl_physicaliam_help.html |pvl_physicaliam|> - Determine the incidence angle modifier based on Snell’s Law
% * <pvl_ashraeiam_help.html |pvl_ashraeiam|> - Determine the incidence angle modifier using the ASHRAE incident angle model
% * <pvl_martinruiziam_help.html |pvl_martinruiziam|> - Determine the incidence angle modifier using the Martin and Ruiz incident angle model
% * <pvl_FSspeccorr_help.html |pvl_FSspeccorr|> - Calculate spectral mismatch modifier based on precipitable water and absolute airmass
% * <pvl_calcparams_desoto_help.html |pvl_calcparams_desoto|> - Create module performance coefficient structure for the single diode model described by DeSoto et al., 2006
% * <pvl_calcparams_CEC_help.html |pvl_calcparams_CEC|> - Create module performance coefficient structure for the single diode model used by the CEC
% * <pvl_calcparams_PVsyst_help.html |pvl_calcparams_PVsyst|> - Create module performance coefficient structure for the single diode model in PVsyst version 6
% * <pvl_singlediode_help.html |pvl_singlediode|> - Solves the single diode equation to obtain a photovoltaic IV curve
% * <pvl_sapm_help.html |pvl_sapm|> - Sandia Array Performance Model to get 5 points on IV curve
% * <pvl_huld_help.html |pvl_huld|> - Calculates DC power using the Huld PV module model
% * <pvl_snlinverter_help.html |pvl_snlinverter|> - Converts DC power and voltage to AC power using Sandia's Grid-Connected PV Inverter model
% * <pvl_adrinverter_help.html |pvl_adrinverter|> - Converts DC power and voltage to AC power using Anton Driesse's Grid-Connected PV Inverter model
% * <pvl_singleaxis_help.html |pvl_singleaxis|> - Determine the rotation angle of a 1 axis tracker, and sun incident angle to tracked surface 
%
%%  Functions for parameter estimation for PV module models
%%
% * <pvl_PVsyst_parameter_estimation_help.html |pvl_PVsyst_parameter_estimation|> - Estimates parameters for the PVsyst version 6 module performance model
% * <pvl_desoto_parameter_estimation_help.html |pvl_desoto_parameter_estimation|> - Estimates parameters for the De Soto single diode module performance model
% * <pvl_huld_parameter_estimation_help.html |pvl_huld_parameter_estimation|> - Estimates parameters for the Huld module performance model
% * <pvl_rectify_IV_curve_help.html |pvl_rectify_IV_curve|> - Ensures that Isc and Voc are included in a IV curve and removes duplicate voltage and current points
% * <pvl_huld_parameter_estimation_help.html |pvl_huld_parameter_estimation|> - Estimates parameters for the Huld module performance model
%
%%  Numerical utilities
%%
% * <pvl_lambertw_help.html |pvl_lambertw|> - Compute values for the Lambert W function W(z)
% * <numdiff_help.html |numdiff|> - Compute numerical derivatives for unequally spaced data
% * <pvl_robustfit_help.html |pvl_robustfit|> - Robust regression for linear models

%% Example Scripts
%%
% * <PVL_TestScript1.html |Example Script 1|> - Example script that
% simulates PV system output for a fixed tilt system using weather data
% from a TMY3 file.
% * <PVL_TestScript2.html |Example Script 2|> - Example script that estimates
% irradiance components, direct normal and diffuse horizontal (DNI and DHI) from global horizonal irradiance
% (GHI)

%%
% Copyright 2015 Sandia National Laboratories


