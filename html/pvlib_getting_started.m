%% PV_LIB Toolbox: Getting Started
% The purpose of the PV Performance Modeling Library (PV_LIB) Toolbox for Matlab 
% is to provide a component based library of open source functions for PV modelers 
% to use to learn about and construct custom models of PV performance.  The initial 
% release of the toolbox includes a complete set of functions that allow a user to 
% convert typical weather (e.g., irradiance, temperature, wind speed, and air pressure) 
% data along with a basic PV system design and simulate the AC power output from a PV system.
%% Standard PV Modeling Steps
% Sandia National Laboratories has defined nine general modeling steps that are usually 
% involved in modeling a PV system’s performance.  These steps assume that the PV system 
% has been defined.  The steps are listed below:

%%
% # *Irradiance and Weather* – Collection of irradiance and weather measurements at standard orientations.
% # *Incident Irradiance* – Transposition and translation of irradiance components to the plane-of-array
% # *Shading and Soiling* -  Reduction in irradiance reaching the module and array
% # *Cell Temperature* – Calculation of cell temperature from irradiance, air temperature, and wind speed.
% # *Module Output* – Determination of the IV curve for the module in the array environment
% # *DC and Mismatch Losses* – Accounting for the effects of mismatch in series connected strings and DC wiring losses.
% # *DC to DC Max Power Point Tracking* - Calculation of the DC voltage required to maximize power.
% # *DC to AC Conversion* – Conversion of the DC power to AC power, including variable efficiency and parasitic losses.
% # *AC Losses* – Accounting for AC losses before the utility meter.

%% What Data is Available?
% The process of modeling the performance of a PV system depends on what kind of weather data 
% is available.  Some typical starting points are listed below.  Example scripts using 
% PV_LIB functions have been developed for these various scenarios.
%%
% * *Standard Weather Data* - (e.g., Typical Meteorological Year)
%% 
% # Direct Normal Irradiance (DNI)
% # Global Horizontal Irradiance (GHI)
% # Diffuse Horizontal Irradiance (DHI)
%
%%
% <PVL_TestScript1.html PVL_TestScript1> provides a complete example of how to 
% model a PV system using standard weather inputs from a TMY3 file.
%%
% * *Global Horizontal Irradiance only*
%%
% <PVL_TestScript2.html PVL_TestScript2> provides an example of how to 
% estimate plane of array beam and diffuse irradiance starting with only
% global horizontal irradiance data.
%%
% Copyright 2015 Sandia National Laboratories
