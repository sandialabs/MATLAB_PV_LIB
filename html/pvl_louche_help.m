%% pvl_louche
% Estimate DNI and DHI from GHI using the Louche model.
%
%% Syntax
% |DNI, DHI, Kt| = pvl_louche(GHI, Z, doy)|
%
%% Description
% The Louche algorithm estimates direct normal irradiance (DNI) from 
% global horizontal irradiance (GHI) using an empirical relationship
% between the direct transmittance and the ratio between GHI and
% extraterrestrial irradiance. DHI is calculated from GHI and DNI as
%
% $$DHI = GHI - DNI \times \cos\left( Z \right) $$
%
% where $Z$ is the true solar zenith angle.
%
%
%% Inputs
% * *|GHI|* - a scalar or vector of global horizontal irradiance in W/m^2.
%     If |GHI| is a vector it must be of the same size as all other vector inputs. 
%     |GHI| must be >=0.
% * *|Z|* - a scalar or vector of true (not refraction-corrected) zenith
%     angles in decimal degrees. If |Z| is a vector it must be of the
%     same size as all other vector inputs. |Z| must be >=0 and <=180.
% * *|doy|* - a scalar or vector of values providing the day of the year. If
%     |doy| is a vector it must be of the same size as all other vector inputs.
%     |doy| must be >= 1 and < 367.
%
%% Outputs
%%
% * *|DNI|* - the modeled direct normal irradiance in W/m^2. 
% * *|DHI|* - the modeled diffuse horizontal irradiance in W/m^2. 
% * *|Kt|* - the ratio of GHI to extraterrestrial irradiance on a horizontal plane.
%
%% Example
% This example compares DNI for August 6 from the TMY3 file for Albuquerque, NM,
% with DNI estimated from GHI using the Louche model.
TMYData = pvl_readtmy3('723650TY.csv');
TimeMatlab = TMYData.DateNumber;
Time = pvl_maketimestruct(TimeMatlab, ones(size(TimeMatlab))*TMYData.SiteTimeZone); 
dayofyear = pvl_date2doy(Time.year, Time.month, Time.day);
DNI = TMYData.DNI; % Read in for comparison with results
DHI = TMYData.DHI; % Read in for comparison with results
GHI = TMYData.GHI;
Location = pvl_makelocationstruct(TMYData.SiteLatitude,TMYData.SiteLongitude,...
TMYData.SiteElevation); %Altitude is optional
PresPa = TMYData.Pressure*100; %Convert pressure from mbar to Pa

% Models that rely on sun position that are run at specific time steps
% (e.g., hourly) run into numerical problems during timesteps when the sun
% straddles the horizon (i.e., spends part of the timestep above and part
% below the horizon).  Since in this example we use TMY day, which is
% hourly and reported at the end of the hour, we adjust the sun position in
% the following two ways.  (1) For hours when the sun is above the horizon
% we adjust sun position so it is for the middle of the hour.  For example,
% sun position reported at 4PM is the sun position at 3:30 PM.  (2) For
% hours where the sun straddles the horizon, we report the position half
% way between the horizon and the position of the sun at the end of the
% hour.  This is appropriate because the GHI data is essentially a sum of the GHI
% measurements made in the past hour.
%
% Run sun position twice:  Once for end of hour positions.  
% A second time for mid hour positions.  
% Adjust sun elevation and GHI for hours when the sun traverses the horizon
% Include these adjusted values in the mid hour positions
%
[~, ~, AppSunEl, ~] = pvl_ephemeris(Time,Location,PresPa,TMYData.DryBulb);
Time.hour = Time.hour-.5; % shift times back 1/2 hour for sun position calculation because of tmy
% timestamps
[~, ~, AppSunEla, ~] = pvl_ephemeris(Time,Location,PresPa,TMYData.DryBulb);

A=diff(sign(AppSunEl)); %identifies hour before sun straddles horizon (2,-2)
ind1 = find(A==2)+1; % AM hour where sun straddles horizon
ind2 = find(A==-2)+1; % PM hour where sun straddles horizon
%AM Adjustment
    AppSunEl(ind1) = AppSunEl(ind1)/2; %change sun elevation to mid way above the horizon 
 %PM Adjustment
    AppSunEl(ind2) = AppSunEl(ind2)/2; %change sun elevation to mid way above the horizon 
AppSunEla(ind1)=  AppSunEl(ind1);
AppSunEla(ind2)=  AppSunEl(ind2);

[DNI_model, DHI_model, ~] = pvl_louche(GHI,90-AppSunEla, dayofyear);

figure
tfilter = and(Time.month == 8,Time.day == 6);
plot(Time.hour(tfilter),DNI_model(tfilter),'-s')
hold all
plot(Time.hour(tfilter),DNI(tfilter),'-o')
xlim([0 25])
legend('DNI (Louche Model)','DNI (TMY3)','Location','NE')
xlabel('Hour of Day')
ylabel('Irradiance (W/m^2)')
title({'Direct Normal Irradiance Comparison';'Albuquerque - Aug 6'},'FontSize',14)
%% References
%%
% * [1] Louche A, Notton G, Poggi P, Simmonnot G. Correlations for direct
% normal and global horizontal irradiation on French Mediterranean site.
% Solar Energy 1991;46:261-6
%   
%

%% See also
% <pvl_disc_help.html |pvl_disc|> , 
% <pvl_dirint_help.html |pvl_dirint|> , 
% <pvl_erbs_help.html |pvl_erbs|> , 
% <pvl_orgill_hollands_help.html |pvl_orgill_hollands|> , 
% <pvl_reindl_1_help.html |pvl_reindl_1|> , 
% <pvl_reindl_2_help.html |pvl_reindl_2|> , 
% <pvl_ephemeris_help.html |pvl_ephemeris|> ,
% <pvl_date2doy_help.html |pvl_date2doy|> ,
% <pvl_alt2pres_help.html |pvl_alt2pres|> 
%%
% Copyright 2014 Sandia National Laboratories
