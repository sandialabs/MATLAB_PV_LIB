%% pvl_Purdue_albedo_model
% Calculate the collection of ground-reflected albedo light on the rear surface of a PV module.
%
%% Syntax
% * |pvl_Purdue_albedo_model(SurfTilt, SurfAz, EtoH, Albedo, DHI, DNI, HExtra, SunZen, SunAz, AM)|
% * |pvl_Purdue_albedo_model(SurfTilt, SurfAz, EtoH, Albedo, DHI, DNI, HExtra, SunZen, SunAz, AM, model)|
%
%% Description
% This code is part of the Purdue Bifacial irradiance model [1] and it can 
% simulate the albedo light intensity on both the front and rear sides of a 
% bifacial solar module. This model takes two types of self-shading losses into
% account: 1) direct blocking of direct beam and circumsolar light by the module onto the ground
% and 2) sky masking of isotropic diffuse light by the module. This model
% employs a view-factor based approach and the detailed methodology is discussed
% in [1].
%
%% Inputs
% * *|SurfTilt|* - a scalar or vector of surface tilt angles in decimal degrees.
%     If |SurfTilt| is a vector it must be of the same size as all other vector
%     inputs. |SurfTilt| must be >=0 and <=180. The tilt angle is defined as
%     degrees from horizontal (e.g. surface facing up = 0, surface facing
%     horizon = 90).
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
%     |DNI| is a vector it must be of the same size as all other vector inputs.
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
%
%% Output
% * *|I_Alb|* - the total ground-reflected albedo irradiance incident to the specified surface.
%     |I_Alb| is a column vector vector with a number of elements equal to the input vector(s).
%
%% References
%   [1] Sun, X., Khan, M. R., Alam, M. A., 2018. Optimization and performance 
%   of bifacial solar modules: A global perspective. Applied Energy 212, pp. 1601-1610.
%   [2] Khan, M. R., Hanna, A., Sun, X., Alam, M. A., 2017. Vertical bifacial solar farms:
%   Physics, design, and global optimization. Applied Energy, 206, 240–248.
%   [3] Duffie, J. A., Beckman, W. A. 2013. Solar Engineering of Thermal Processes (4th Editio). 
%   Wiley.
%
%% See also 
% <pvl_perez_help.html |pvl_perez|> ,
% <pvl_Purdue_bifacial_irradiance_help.html |pvl_Purdue_bifacial_irradiance|>
%
%% Notes
% |pvl_Purdue_albedo_model| contributed by Xingshu Sun of Purdue
% University, 2018.


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

%% Describe south-facing, latitude tilt module
SurfTilt_Front = 30;
SurfTilt_Rear = 180 - SurfTilt_Front;
SurfAz_Front = 180;
SurfAz_Rear = 180 + SurfAz_Front;
Albedo = 0.25; % 25% albedo coefficient for concrete/vegetation

Module_length = 1.4; % meters, dimension of module not parallel to the ground
Module_elev = [1 0]; % meters, distance from the ground to the bottom of the module

EtoH = Module_elev./Module_length; % Elevated (EtoH > 1) and flush with the ground (EtoH = 0) 

% run calculation
for i = 1:length(EtoH);
    
    I_Alb_Front(:,i) = pvl_Purdue_albedo_model(SurfTilt_Front, SurfAz_Front, EtoH(i), Albedo, ...
        TMYData.DHI(tfilter), TMYData.DNI(tfilter), HExtra, SunZen, SunAz, AM);
    I_Alb_Rear(:,i) = pvl_Purdue_albedo_model(SurfTilt_Rear, SurfAz_Rear, EtoH(i), Albedo, ...
        TMYData.DHI(tfilter), TMYData.DNI(tfilter), HExtra, SunZen, SunAz, AM);
end


%%

figure
hold all
s = {'r-','b-','ro','bo'};
for i=1:length(EtoH);
    plot(Time.hour, I_Alb_Rear(:,i), s{2*i})
    plot(Time.hour, I_Alb_Front(:,i), s{2*i-1})
end

xlim([-5 22])
xlabel('Hour of Day')
ylabel('Albedo Irradiance (W/m^2)')
legend('Rear (elevated)', 'Front (elevated)', 'Rear (Ground Mounted)', 'Front (Ground Mounted)', ...
    'Location', 'NorthWest')
title({'Albedo Irradiance on a Bifacial Module';'Albuquerque - Aug 2'},'FontSize',14)
