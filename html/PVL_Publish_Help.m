% Publish all help files
addpath(genpath('S:\MATLAB_Solar_Functions\PV_LIB Version 1_3 Development'))

%% Time and Location Utilities
publish('pvl_date2doy_help.m');             %1
publish('pvl_doy2date_help.m');             %2
publish('pvl_leapyear_help.m');             %3
publish('pvl_matlabtime2excel_help.m');     %4
publish('pvl_exceltime2matlab_help.m');     %5
publish('pvl_maketimestruct_help.m');       %6
publish('pvl_makelocationstruct_help.m');   %7
%% Irradiance and Atmospheric Functions
publish('pvl_readtmy2_help.m');             %8
publish('pvl_readtmy3_help.m');             %9
publish('pvl_getISDdata_help.m');
publish('pvl_readISH_help.m');
publish('pvl_ephemeris_help.m');            %10
publish('pvl_spa_help.m');                  %11
publish('pvl_extraradiation_help.m');       %12
publish('pvl_pres2alt_help.m');             %13
publish('pvl_alt2pres_help.m');             %14
publish('pvl_relativeairmass_help.m');      %15
publish('pvl_absoluteairmass_help.m');      %16
publish('pvl_disc_help.m');                 %17
publish('pvl_dirint_help.m');               %18
publish('pvl_louche_help.m');               %19
publish('pvl_orgill_hollands_help.m');      %20
publish('pvl_erbs_help.m');                 %21
publish('pvl_reindl_1_help.m');             %22
publish('pvl_reindl_2_help.m');             %23
publish('pvl_clearsky_haurwitz_help.m');    %24
publish('pvl_clearsky_ineichen_help.m');    %25
publish('pvl_calcPwat_help.m');

%% Irradiance Translation Functions
publish('pvl_grounddiffuse_help.m');        %26
publish('pvl_isotropicsky_help.m');         %27
publish('pvl_reindl1990_help.m');           %28
publish('pvl_perez_help.m');                %29
publish('pvl_kingdiffuse_help.m');          %30
publish('pvl_klucher1979_help.m');          %31
publish('pvl_haydavies1980_help.m');        %32
publish('pvl_getaoi_help.m');               %33

%% Photovoltaic System Functions
publish('pvl_sapmmoduledb_help.m');         %37
%publish('pvl_snlinverterdb_help.m');        %42
publish('pvl_SAMLibraryReader_CECModules_help.m');
publish('pvl_SAMLibraryReader_SNLInverters_help.m');
publish('pvl_sapmcelltemp_help.m');         %38

publish('pvl_physicaliam_help.m');          %34
publish('pvl_ashraeiam_help.m');            %35
publish('pvl_martinruiziam_help.m');        %36

publish('pvl_FSspeccorr_help.m');

publish('pvl_calcparams_desoto_help.m');    %40
publish('pvl_calcparams_CEC_help.m');
publish('pvl_calcparams_PVsyst_help.m')
publish('pvl_singlediode_help.m');          %41

publish('pvl_sapm_help.m');                 %39
publish('pvl_huld_help.m');
publish('pvl_snlinverter_help.m');          %43
publish('pvl_adrinverter_help.m');          %44

publish('pvl_singleaxis_help.m');           %44

publish('pvl_PVsyst_parameter_estimation_help.m');
publish('pvl_desoto_parameter_estimation_help.m');
publish('pvl_rectify_IV_curve_help.m');

%% Numerical utilities
publish('pvl_lambertw_help.m');
publish('numdiff_help.m');

%% Example Scripts
publish('PVL_TestScript1.m');               
publish('PVL_TestScript2.m');               

%% Other Help
publish('pvlib_features.m');
publish('pvlib_functions_by_cat.m');
publish('pvlib_getting_started.m');
publish('pvlib_product_page.m');
publish('pvlib_release_notes.m');
publish('pvlib_system_requirements.m');
publish('pvlib_user_guide.m');



