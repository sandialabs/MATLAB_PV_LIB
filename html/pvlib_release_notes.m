%% PV_LIB Toolbox version 1.3 Release Notes
% This toolbox implements functions that enable simulation of the
% performance of photovoltaic (PV) energy systems.
%
% The PV_LIB Toolbox requires Matlab software.
% 
%% References
% Additional information and documentation is available on the PV
% Performance Modeling Collaborative website (<http://pvpmc.org>).
%
%
%% Bug Reporting
% Report bugs and problems to Joshua Stein (jsstein@sandia.gov).
%
%% Credits for Non-Sandia contributions
% * *Mitchell Lee and Alex Panchula, FirstSolar:* Contributed the spectral
% mismatch modifier model in functions |pvl_FSspeccorr| and |pvl_calcPwat|.
% * *Rob Andrews, Queen's University:* Multiple bug finds and fixes in
% PV_LIB version 1.0 functions |pvl_perez|, |pvl_haydavies1980|,
% |pvl_klucher1979|, and |pvl_reindl1990|.
% * *Martin Herrerias Azcue:* Found incorrect documentation notes for
% |pvl_calcparams_desoto|.
% * *John McKeen, DOW Solar Solutions:* Found incorrect error checking in
% |pvl_calcparams_desoto|.
% * *Mark Campanelli, NREL:* Found bug in |pvl_clearsky_ineichen|.
%
%% Version History 
%
%
%% *Version 1.31:* January-2016
%
% # Includes |est_single_diode_param|, |Schumaker_QSpline|, 
% |calc_theta_phi_exact|, |update_Io_known_n|, and |update_Rsh_fixed_pt|
% which were left out of the Version 1.3 release.
% # Removed dependency on the Statistics toolbox from
% |est_single_diode_param|.
%
%% *Version 1.3:* December-2015
%
% *New Functions*
%
% |pvl_PVsyst_parameter_estimation| - Estimates parameters for the PVsyst
% version 6 module performance model.
%
% |pvl_desoto_parameter_estimation| - Estimates parameters for the De Soto
% single diode module performance model.
%
% |pvl_rectify_IV_curve| - Ensures that Isc and Voc are included
% in a IV curve and removes duplicate voltage and current points.
%
% |pvl_calcparams_CEC| - Create module performance model coefficients
% for the single diode model used by the CEC.
%
% |pvl_calcparams_PVsyst| - Create module performance model coefficient
% for the single diode model in PVsyst version 6.
%
% |pvl_huld| - Calculates DC power using the Huld PV module model.
%
% |pvl_calcPwat| - Calculate precipitable water (cm) from ambient
% air temperature (C) and relatively humidity (%).
%
% |pvl_FSspeccorr| - Calculate spectral mismatch modifier based on
% precipitable water and absolute airmass.
%
% |pvl_getISDdata| - Fetch data from NOAA ISD at ftp.ncdc.noaa.gov.
%
% |pvl_readISH| - Read data fetched from ftp.ncdc.noaa.gov into a table.
%
% |pvl_lambertw| - Compute values for the Lambert W function W(z).
%
% |numdiff| - Compute numerical derivatives for unequally spaced data.
%
% *Significant Changes*
%
% # Corrected |pvl_erbs|; the coefficient of K_t^2 changed from 4.3888 to
% 4.388
% # Changed |pvl_calcparams_desoto| to require effective irradiance as its
% input, rather than absorbed irradiance (i.e., broadband irradiance
% reaching a module's cells). |pvl_calcparams_desoto| no longer computes and
% applies a spectral mismatch modifier.
% # Changed |pvl_calcparams_desoto| to remove the requirement that all inputs
% be vectors of the same length.  Inputs |S| and |Tcell| can be vectors of the
% same length; all other inputs should be scalars.
% # Replaced |wapr_vec| with |pvl_lambertw|. |wapr_vec| is obsolete and removed.
% # Changed |pvl_haydavies1980| and |pvl_perez| to accept |NaN| in inputs.
% # Added functions |pvl_FSspeccorr| and |pvl_calcPwat| to implement
% FirstSolar’s spectral mismatch modifier model.
% # Removed |pvl_snlinverterdb| (made obsolete by
% |pvl_SAMLibraryReader_SNLInverters| and |SandiaInverterDatabaseSAM2014.1.14.mat|).
% # Changed |pvl_erbs|, |pvl_louche|, |pvl_orgill_hollands|, |pvl_reindl_1|,
% and |pvl_reindl_2| to output DNI, DHI and Kt, and to clarify each model’s description.
% # Corrected description of |pvl_kingdiffuse| to clarify that the King 
% diffuse model estimates total plane-of-array diffuse (sky diffuse and ground reflected irradiance).
% # Included databases of CEC module model coefficients (|CECModuleDatabaseSAM2015.6.30.mat|) and Sandia inverter
% model coefficients (|CECInverterDatabaseSAM2015.6.30.mat\) from SAM2015.6.30. 
%
%% *Version 1.2:* December-2014
%
% *New Functions*
%
% |pvl_louche| - Estimate Direct Normal Irradiance from Global Horizontal
% Irradiance using the Louche model.
%
% |pvl_erbs| - Estimate Direct Normal Irradiance from Global Horizontal
% Irradiance using the Erbs model.
%
% |pvl_orgill_hollands| - Estimate Direct Normal Irradiance from Global
% Horizontal Irradiance using the Orgill and Hollands model.
%
% |pvl_reindl_1| - Estimate Direct Normal Irradiance from Global Horizontal
% Irradiance using the Reindl 1 model.
%
% |pvl_reindl_2| - Estimate Direct Normal Irradiance from Global Horizontal
% Irradiance using the Reindl 2 model.
%
% |pvl_martinruiziam| - Determine the incidence angle modifier using the
% Martin and Ruiz incident angle model.
%
% |pvl_adrinverter| - Converts DC power and voltage to AC power using Anton
% Driesse's Grid-Connected PV Inverter model.
%
% *Significant Changes*
%
% * Corrected |pvl_klucher1979|. Statement |GHI(GHI<DHI) = DHI| did not
% work for vectors, changed to |GHI(GHI<DHI) = DHI(GHI<DHI)|.
% * Changed |pvl_calcparams_desoto.m| as follows:
%
% # Removed line |M = max(M, 0)|.
% # Removed line |S(S==0) = 1E-10|. 
% # Added line |IL(isnan(M) | M<0 | S <=0) = 0|. This sets IL to 0 when 
% irradiance is <0, airmass modifier is <0 (most likely due to evaluating 
% the polynomial at very high airmass or airmass modifier is |NaN|, which
% is returned by |pvl_relativeairmass| for sun zenith angles > 90 degrees.
% # Added line |I0(IL==0) = 0|. According the circuit diagram, IL is the only source, and therefore, there can be no reverse saturation current (I0) when IL is 0
% # Added line |Rsh(S <= 0) = inf|. This is due to the fact that Rsh is determined by dividing by S, negative values would give a negative Rsh, and at S=0 Rsh is undefined
%
% * Modified |pvl_singlediode.m| 
%
% # Added a filter |u = IL > 0|. Thus we only compute IV points and IV
% curves when there is a photocurrent (IL).
% # Pre-allocate memory to Imax, negPmp, Vmax, Ix, Ixx, Voc, and Isc.
% This is necessary in order to use the filter.
% # Changed the Isc and Voc generation lines to only generate Isc and Voc under conditions which satisfy filter u.
% # Changed Imax and negPmp finding line to only find Imax and negPmp
% when conditions satisfy filter u.
% # Changed Vmax, Ixx, and Ix generation lines to only generate the
% values under conditions which satisfy filter u.
% # Pre-allocate memory for Result.I and Result.V (necessary to implement filter u).
% # Changed |Result.I| finding line in order to only find I-V curve
% currents under conditions which satisfy filter u.
% # Added line |Result.I(:,end) = 0| in order to ensure that the current at Voc is exactly 0 (prevents numerical error which may result in a negative current at Voc).
%
% * Corrected errors in |wapr_vec.m|.  From line 134 through line 154 (inclusive), the Lambert W function is evaluated piecewise over the domain (-exp(-1), inf) using four approximations.  The domain is partitioned into four disjoint sets and each approximation is applied to one (and only one) set.  Filters were implemented incorrectly which allowed later approximations to overwrite values obtained from early approximations.
% * Corrected |pvl_calcparams_desoto.m| documentation. The documentation used to specify that alpha_isc should be in units of 1/C, it now specifies that alpha_isc should be in units of A/C or A/K. Documentation was changed rather than code due to the fact that the SAM CEC module database lists alpha_isc in units of A/C. Thanks to *Martin Herrerias Azcue* for the bug find.
% * Corrected |pvl_calcparams_desoto.m|. Prior input checking for |EgRef| was
% |p.addRequired('EgRef', @(x) (isnumeric(x) & isvector(x) & x>0))|; this has
% been changed to |p.addRequired('EgRef', @(x) (isnumeric(x) & isvector(x) & all(x>0)))|;
% in order to ensure that all |EgRef| values are positive. Thanks to 
% *John McKeen* at DOW Solar Solutions for the bug find.
% * Corrected |pvl_clearsky_ineichen| to use system-dependent file separators
% when using the default Linke turbidity index. Prior code was
% |load('Required Data\LinkeTurbidities.mat');| this has been changed to 
% |load(['Required Data' filesep 'LinkeTurbidities.mat']);| in order to allow for system-specific file separators.
% Thanks to *Mark Campanelli* at NREL for finding this bug.
% * Modified |pvl_singlediode| to not require |pvl_fminbnd| in order to find
% the maximum power point. Prior to this change, the maximum power point 
% was found by finding the maximum value of the power-current curve using 
% a modified version of MATLAB’s |fminbnd| function. |pvl_singlediode| now 
% finds where dP/dV as a function of current is equal to 0. Uses bisection 
% techniques. New subfunctions |calc_phi_exact|, |calc_Imp_bisect|, |g|, and |calc_Pmp_bisect| were added.
% * Modified |pvl_getaoi| to avoid complex values. Changed the line  
% |AOI = acosd(cosd(SunZen).*cosd(SurfTilt)+sind(SurfTilt).*sind(SunZen).*cosd(SunAz-SurfAz));| 
% to |AOI = acosd(max(min(cosd(SunZen).*cosd(SurfTilt)+sind(SurfTilt).*sind(SunZen).*cosd(SunAz-SurfAz),
% 1),-1));|. Roundoff errors previously could cause the argument of the arcos function to be greater than 1 or less than -1 (resulting in a complex output). The min and max functions prevent such an occurrence.
% * Modified |pvl_getaoi| to avoid complex values. Changed the line  
% |AOI = acosd(max(min(cosd(SunZen).*cosd(SurfTilt)+sind(SurfTilt).*sind(SunZen).*cosd(SunAz-SurfAz), 1),-1));|
% to |temp = cosd(SunZen).*cosd(SurfTilt)+sind(SurfTilt).*sind(SunZen).*cosd(SunAz-SurfAz);| 
% |temp(temp>1) = 1; temp(temp<-1) = -1;| 
% |AOI = acosd(temp);|.
% Roundoff errors previously could cause the argument of the arcos function to be greater than 1 or less than -1 (resulting in a complex output). The min and max functions did not prevent all occurrences.
%
%% *Version 1.1:* December-2012
% 
% *New functions*
%
% |pvl_clearsky_haurwitz| - Clear sky GHI model.
%
% |pvl_clearsky_ineichen| - Clear sky irradiance (GHI, DNI, and DHI) model.
%
% |pvl_physicaliam| - Incident angle modifer model based on Snell’s Law.
%
% |pvl_ashraeiam| - Incident angle modifier model from ASHRAE.
%
% |pvl_calcparams_desoto| – Calculate module performance model coefficients
% for the De soto single diode model.
%
% |pvl_singlediode| – Solves the single-diode equation to obtain a
% photovoltaic IV curve.
%
% |pvl_SAMLibraryReader_CECModules| - Open a System Advisor Model (SAM) CEC
% module library.
%
% |pvl_SAMLibraryReader_SNLInverters| - Open a System Advisor Model (SAM)
% Sandia inverter library.
%
% *Significant Changes*
%
% * Fixed numerous text typos in documentation and help files.
% * Made numerical fix to |pvl_spa|.
% * Fixed angle of incidence calculation in |pvl_perez|,
% |pvl_haydavies1980|, |pvl_klucher 1979|, and |pvl_reindl1990|.
% * Fixed numerical errors in final calculation of |pvl_perez|.
% * Fixed |pvl_perez| to accept scalar input values with vector input values.
%
%% *Version 1.0:* June-2012  Initial Release
%
% Copyright 2015 Sandia National Laboratories