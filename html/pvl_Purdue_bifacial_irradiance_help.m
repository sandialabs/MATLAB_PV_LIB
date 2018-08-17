%% pvl_Purdue_bifacial_irradiance
% Calculate the irradiance on the front and rear sides of a bifacial solar module.
%
%% Syntax
% * *|pvl_Purdue_bifacial_irradiance(SurfTilt, SurfAz, EtoH, DHI, DNI, Albedo, HExtra, SunZen, SunAz, AM)|
% * *|pvl_Purdue_bifacial_irradiance(SurfTilt, SurfAz, EtoH, DHI, DNI, Albedo, HExtra, SunZen, SunAz, AM, model)|
% * *|pvl_Purdue_bifacial_irradiance(SurfTilt, SurfAz, EtoH, DHI, DNI, Albedo, HExtra, SunZen, SunAz, AM, model, Rloss)|
% * *|pvl_Purdue_bifacial_irradiance(SurfTilt, SurfAz, EtoH, DHI, DNI, Albedo, HExtra, SunZen, SunAz, AM, model, Rloss, Rloss_Para)|
%
%% Description
%   The Purdue Bifacial irradiance model [1] simulates the total irradiance including 
%   direct, diffuse, and albedo light, on both the front and rear sides of
%   a bifacial solar module. This model applies an analytical view-factor 
%   based approach to explicitly account for the self-shading losses of 
%   albedo light due to (1) direct blocking of direct and circumsolar diffuse
%   light and 2) sky masking of isotropic diffuse light onto the ground. 
%   This model also incorporates an optional reflection-loss model [4]. 
%   This model has been validated against data spanning from Africa, Europe,
%   Asia, and North America (please refer [1] for more detail).
%
%% Inputs
% * *|SurfTilt|* - a scalar or vector of surface tilt angles in decimal degrees.
%     The tilt angle is defined as degrees from horizontal (e.g. surface facing
%     up = 0, surface facing horizon = 90). |SurfTilt| must be >=0 and <=180.
%     If |SurfTilt| is a vector it must be of the same size as all other vector inputs.  
% * *|SurfAz|* - a scalar or vector of surface azimuth angles in decimal degrees.
%     If |SurfAz| is a vector it must be of the same size as all other vector
%     inputs. |SurfAz| must be >=0 and <=360. The Azimuth convention is defined
%     as degrees east of north (e.g. North = 0, East = 90, West = 270).
% * *|EtoH|* - a scalar or vector of the ratio of module elevation(E) to module height(H).
%     Module height is the module dimension not parallel to the ground.
%     If |EtoH| is a vector it must be of the same size as all other vector
%     inputs. |EtoH| must be >=0.
% * *|Albedo|* - a scalar or vector of groud albedo coefficient.
%     If |Albedo| is a vector it must be of the same size as all other vector
%     inputs. |Albedo| must be >=0 and <=1.
% * *|DHI|* - a scalar or vector of diffuse horizontal irradiance in W/m^2. 
%     If |DHI| is a vector it must be of the same size as all other vector inputs.
%     |DHI| must be >=0.
% * *|DNI|* - a scalar or vector of direct normal irradiance in W/m^2. If
% * *|DNI| is a vector it must be of the same size as all other vector inputs.
%     |DNI| must be >=0.
% * *|HExtra|* - a scalar or vector of extraterrestrial normal irradiance in
%     W/m^2. If |HExtra| is a vector it must be of the same size as
%     all other vector inputs. |HExtra| must be >=0.
% * *|SunZen|* - a scalar or vector of apparent (refraction-corrected) zenith
%     angles in decimal degrees. If |SunZen| is a vector it must be of the
%     same size as all other vector inputs. |SunZen| must be >=0 and <=180.
% * *|SunAz|* - a scalar or vector of sun azimuth angles in decimal degrees.
%     If |SunAz| is a vector it must be of the same size as all other vector
%     inputs. |SunAz| must be >=0 and <=360. The Azimuth convention is defined
%     as degrees east of north (e.g. North = 0, East = 90, West = 270).
% * *|AM|* - a scalar or vector of relative (not pressure-corrected) airmass 
%     values. If |AM| is a vector it must be of the same size as all other 
%     vector inputs. |AM| must be >=0.
% * *|model|* - a character string which selects the desired set of Perez
%     coefficients. If model is not provided as an input, the default,
%     '1990' will be used.
%     All possible model selections are: 
%       '1990', 'allsitescomposite1990' (same as '1990'),
%       'allsitescomposite1988', 'sandiacomposite1988',
%       'usacomposite1988', 'france1988', 'phoenix1988',
%       'elmonte1988', 'osage1988', 'albuquerque1988',
%       'capecanaveral1988', or 'albany1988'
% * *|Rloss|* - a character string which determines the inclusion of reflection 
%     loss model. By default, 'Yes' will be used. If 'No' is input, reflection 
%     loss will be neglected.
% * *|Rloss_Para|* - a three-element vector represents the parameters 
%     (in order, ar, c1, and c2) in the reflection models in Ref. [4].
%     By default a parameter set for a glass-faced silicon solar module, 
%     [ar = 0.16, cl = 4/(3*pi), c2 = -0.074], will be used.  
%
%% Output
% * *|Front_Irradiance|* - the total irradiance including direct beam, diffuse,
%     and albedo light on the front side. |Front_Irradiance| is a column vector
%     of the same size as the input vector(s).
% * *|Rear_Irradiance|* - the total irradiance includig direct beam, diffuse,
%     and albedo light on the rear side.  |Rear_Irradiance| is a column vector
%     of the same size as the input vector(s).
%
%% References
%   [1] Sun, X., Khan, M. R., Alam, M. A., 2018. Optimization and performance 
%   of bifacial solar modules: A global perspective. Applied Energy 212, pp. 1601-1610.
%   [2] Khan, M. R., Hanna, A., Sun, X., Alam, M. A., 2017. Vertical bifacial solar farms:
%   Physics, design, and global optimization. Applied Energy, 206, 240–248.
%   [3] Duffie, J. A., Beckman, W. A. 2013. Solar Engineering of Thermal Processes (4th Editio). 
%   Wiley.
%   [4] Martín, N., Ruiz, J. M. 2005. Annual angular reflection losses in
%   PV modules. Progress in Photovoltaics: Research and Applications, 13(1), 75–84.
%
%% See also 
% <pvl_perez_help.html |pvl_perez|> ,
% <pvl_Purdue_albedo_model_help.html |pvl_Purdue_albedo_model|> ,
% <pvl_iam_martinruiz_components_help.html |pvl_iam_martinruiz_components|>
%
%% Notes
% |pvl_Purdue_bifacial_irradiance| contributed by Xingshu Sun of Purdue
% University, 2018.

%% Example
clearvars
close all

%Load weather data
TMYData = pvl_readtmy3('723650TY.csv');
TimeMatlab = TMYData.DateNumber;
dv = datevec(TimeMatlab);

tfilter = and(dv(:,2) == 8, dv(:,3) == 2); % Select August 2

Time = pvl_maketimestruct(TimeMatlab(tfilter), ones(size(TimeMatlab(tfilter)))*TMYData.SiteTimeZone);

% Sun position calculations
HExtra = pvl_extraradiation(pvl_date2doy(Time.year,Time.month,Time.day));
Location = pvl_makelocationstruct(TMYData.SiteLatitude,TMYData.SiteLongitude,TMYData.SiteElevation);
PresPa = TMYData.Pressure(tfilter)*100; %Convert pressure from mbar to Pa
[SunAz, SunEl, AppSunEl, SolarTime] = pvl_ephemeris(Time,Location,PresPa,TMYData.DryBulb(tfilter));
SunZen = 90 - AppSunEl;
AM = pvl_relativeairmass(SunZen);
AM(isnan(AM)) = 20;

%% Describe vertical east-west module
Module_length = 1.4; % meters, dimension of module not parallel to the ground
Module_elev = 0.5; % meters, distance from the ground to the bottom of the module
SurfTilt = 90; %Vertical installation
SurfAz = 90; %East-west facing
EtoH = Module_elev / Module_length; %Ground-mounted
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
title({'Irradiance on a east-facing vertical bifacial module';'Albuquerque - Aug 2'},'FontSize',14)


%% Describe south-facing latitude tilt module mounted 0.5m from the ground
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

                  