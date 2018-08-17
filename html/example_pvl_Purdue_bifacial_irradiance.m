%% example_pvl_Purdue_albedo_model
%
% Example calculation of front and rear irradiance on a bifacial module.

clc
clearvars
close all

%% Load weather data
TMYData = pvl_readtmy3('723650TY.csv');
TimeMatlab = TMYData.DateNumber;
dv = datevec(TimeMatlab);

tfilter = and(dv(:,2) == 8, dv(:,3) == 2); % Select August 2

Time = pvl_maketimestruct(TimeMatlab(tfilter), ones(size(TimeMatlab(tfilter)))*TMYData.SiteTimeZone);

%% Sun position calculations
HExtra = pvl_extraradiation(pvl_date2doy(Time.year,Time.month,Time.day));
Location = pvl_makelocationstruct(TMYData.SiteLatitude,TMYData.SiteLongitude,TMYData.SiteElevation);
PresPa = TMYData.Pressure(tfilter)*100; %Convert pressure from mbar to Pa
[SunAz, SunEl, AppSunEl, SolarTime] = pvl_ephemeris(Time,Location,PresPa,TMYData.DryBulb(tfilter));
SunZen = 90 - AppSunEl;
AM = pvl_relativeairmass(SunZen);
AM(isnan(AM)) = 20;

%% Describe system
SurfTilt = 90; %Vertical installation
SurfAz = 90; %East-west facing
EtoH = 0; %Ground-mounted
Albedo = 0.25; % 25% albedo coefficient for concrete/vegetation
Rloss_Para = [0.16 4/(3*pi) -0.074]; %Default parameters to calculate the reflection loss between the air/glass interface

%% Run the irradiance calculations
[I_Front_R,I_Rear_R] = pvl_Purdue_bifacial_irradiance(SurfTilt, SurfAz, EtoH, Albedo,...
      TMYData.DHI(tfilter), TMYData.DNI(tfilter), HExtra, SunZen, SunAz, AM,'1990','Yes',Rloss_Para); %Consider reflection loss
  
[I_Front_NR,I_Rear_NR] = pvl_Purdue_bifacial_irradiance(SurfTilt, SurfAz, EtoH, Albedo,...
      TMYData.DHI(tfilter), TMYData.DNI(tfilter), HExtra, SunZen, SunAz, AM,'1990','No'); %Do NOT consider reflection loss

%% Plot the results

figure
hold all
s = {'r-','b-','k-','r--','b--','k--'};
plot(Time.hour,I_Front_NR,s{4})
plot(Time.hour,I_Rear_NR,s{5})
plot(Time.hour,I_Front_NR+I_Rear_NR,s{6})
plot(Time.hour,I_Front_R,s{1})
plot(Time.hour,I_Rear_R,s{2})
plot(Time.hour,I_Front_R+I_Rear_R,s{3})
xlim([-5 22])
xlabel('Hour of Day')
ylabel('Total Irradiance (W/m^2)')
legend('Front','Rear','Total', ...
    'Front red. by refl.','Rear red. by refl.','Total red. by refl.', ...
    'Location', 'NorthWest')
title({'Irradiance on a east-facing vertical bifacial module';'Albuquerque - Aug 2'},'FontSize',14)


%% Describe system
Module_length = 1.4; % meters, dimension of module not parallel to the ground
Module_elev = 0.5; % meters, distance from the ground to the bottom of the module
SurfTilt = 35; % Latitude tilt installation
SurfAz = 180; % South facing
EtoH = Module_elev / Module_length; 

Albedo = 0.25; % 25% albedo coefficient for concrete/vegetation
Rloss_Para = [0.16 4/(3*pi) -0.074]; %Default parameters to calculate the reflection loss between the air/glass interface

%Running the calculation
[I_Front_R,I_Rear_R] = pvl_Purdue_bifacial_irradiance(SurfTilt, SurfAz, EtoH, Albedo,...
      TMYData.DHI(tfilter), TMYData.DNI(tfilter), HExtra, SunZen, SunAz, AM,'1990','Yes',Rloss_Para); %Consider reflection loss
  
[I_Front_NR,I_Rear_NR] = pvl_Purdue_bifacial_irradiance(SurfTilt, SurfAz, EtoH, Albedo,...
      TMYData.DHI(tfilter), TMYData.DNI(tfilter), HExtra, SunZen, SunAz, AM,'1990','No'); %Do NOT consider reflection loss

%% Plot the results

figure
hold all
s = {'r-','b-','k-','r--','b--','k--'};
plot(Time.hour,I_Front_NR,s{4})
plot(Time.hour,I_Rear_NR,s{5})
plot(Time.hour,I_Front_NR+I_Rear_NR,s{6})
plot(Time.hour,I_Front_R,s{1})
plot(Time.hour,I_Rear_R,s{2})
plot(Time.hour,I_Front_R+I_Rear_R,s{3})
xlim([-5 22])
xlabel('Hour of Day')
ylabel('Total Irradiance (W/m^2)')
legend('Front','Rear','Total', ...
    'Front red. by refl.','Rear red. by refl.','Total red. by refl.', ...
    'Location', 'NorthWest')
title({'Irradiance on a south-facing latitude tilt bifacial module';'Albuquerque - Aug 2'},'FontSize',14)

%%
% Copyright 2018 Sandia National Laboratories
