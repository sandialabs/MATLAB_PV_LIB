function SR = pvl_hsu_soiling(Time, Rain, RainThresh, Tilt, PM2_5, PM10, varargin)
% PVL_HSU_SOILING Calculates soiling rate over time given particulate and rain data
%
% Syntax
%   [SR]=pvl_hsu_soiling(Time, Rain, RainThresh, Tilt, PM2_5, PM10)
%   [SR]=pvl_hsu_soiling(Time, Rain, RainThresh, Tilt, PM2_5, PM10, ModelType)
%   [SR]=pvl_hsu_soiling(Time, Rain, RainThresh, Tilt, PM2_5, PM10, ModelType, RainAccPeriod)
%   [SR]=pvl_hsu_soiling(Time, Rain, RainThresh, Tilt, PM2_5, PM10, ModelType, RainAccPeriod, LUC, WindSpeed, Temperature)
%
% Description  
%   [SR]=pvl_hsu_soiling(Time, Rain, RainThresh, Tilt, PM2_5, PM10)
%      Uses the time, rainfall amounts, tilt angle, particulate matter
%      concentrations, and rain cleaning threshold to predict soiling rate
%      using the fixed settling model. Uses hourly rain accumulation to
%      determine cleaning.
%   [SR]=pvl_hsu_soiling(Time, Rain, RainThresh, Tilt, PM2_5, PM10, ModelType)
%      Uses the time, rainfall amounts, tilt angle, particulate matter
%      concentrations, and rain cleaning threshold to predict soiling rate
%      using the model type selected by ModelType. Uses hourly rain
%      accumulation to determine cleaning.
%   [SR]=pvl_hsu_soiling(Time, Rain, RainThresh, Tilt, PM2_5, PM10, ModelType, RainAccPeriod)
%      Uses the time, rainfall amounts, tilt angle, particulate matter
%      concentrations, and rain cleaning threshold to predict soiling rate
%      using the model type selected by ModelType. The user specifies the
%      rain accumulation period (in hours) over which to determine whether
%      sufficient rain has fallen to clean the module.
%   [SR]=pvl_hsu_soiling(Time, Rain, RainThresh, Tilt, PM2_5, PM10, ModelType, RainAccPeriod, LUC, WindSpeed, Temperature)
%      Uses the parameters of the fixed models along with Land Use Category
%      (LUC), wind speed, and ambient temperature to predict soiling ratio
%      using the variable deposition velocity model.
%
% Input Parameters:
%   Time is a struct with elements as described in pvl_maketimestruct. 
%     Time values for the soiling function do not need to be 
%     regularly spaced, although large gaps in timing are discouraged.
%
%   Rain is a vector of rainfall values of the same length as the fields of
%     Time, where each element of Rain correlates with the element of Time.
%     Rainfall values should be in mm of rainfall. Programmatically, rain
%     is accumulated over a given time period, and cleaning is applied
%     immediately after a time period where the cleaning threshold is
%     reached.
%
%   RainThresh is a scalar for the amount of rain, in mm, in a an accumulation
%     period needed to clean the PV modules. In periods where the 
%     accumulated rain meets or exceeds RainThresh, the panels are assumed 
%     to be cleaned immediately after the accumulation period 
%     [1] suggests a good RainThresh could be 1mm, but the time period is
%     not specified. Rain accumulation period length can be adjusted in the
%     optional input RainAccPeriod.
%
%   Tilt is a scalar or vector for the tilt of the PV panels. Changing tilt
%     angles (e.g. in tracking cases) can be accomodated, and tilt angles
%     are correlated with the entries in Time.
%
%   PM2_5 is the concentration of particulate matter with diamter less than
%     2.5 microns, in g/m^3. Historical US concentration data sets may be
%     found at https://aws.epa.gov/aqsweb/airdata/download_files.html. With
%     file descriptions (for hourly data sets) at 
%     https://aqs.epa.gov/aqsweb/airdata/FileFormats.html#_hourly_data_files
%   
%   PM10 is the concentration of particulate matter with diamter less than
%     10 microns, in g/m^3. Historical US PM data sets found as described
%     above.
%
%   ModelType is an optional input to the function to determine the 
%     the model type to be used in the soiling model. A
%     value of "1" indicates that the Variable Deposition Velocity model
%     shall be used, a value of "2" indicates that the Fixed Settling
%     Velocity model shall be used, and a value of "3" indicates that the
%     Fixed Deposition Velocity model shall be used. [1] indicates that the
%     Fixed Settling Velocity model performs best under a wide range of
%     conditions, and thus "2" is the default ModelType if none is entered. 
%     Validation efforts by Sandia National Laboratories
%     confirms these findings. If an incorrect ModelType is provided, the
%     Fixed Settling Velocity will be used (with a warning).
%
%   RainAccPeriod is an optional input that specifies the period, in hours
%     over which to accumulate rainfall totals before checking against the
%     rain cleaning threshold. For example, if the rain threshold is
%     0.5 mm per hour, then RainThresh should be 0.5 and RainAccPeriod
%     should be 1. If the threshold is 1 mm per hour, then the values
%     should be 1 and 1, respectively. The minimum RainAccPeriod is 1
%     hour. The default value is 1, indicating hourly rain accumulation.
%     Accumulation periods exceeding 24 (daily accumulation) are not
%     recommended.
%   
%   LUC is an optional input to the function, but it is required for the
%     Variable Deposition Model. LUC is the Land Use Category as specified
%     in Table 19.2 of [2]. LUC must be a numeric scalar with value 1, 4,
%     6, 8, or 10, corresponding to land with evergreen trees, deciduous
%     trees, grass, desert, or shrubs with interrupted woodlands. The LUC
%     changes the mass deposition by Brownian diffusion. If omitted, the
%     default value of 8 (desert) is used. Note that the value of gamma is
%     the only thing changed by the LUC. LUC 1 and 4 use gamma = 0.54,
%     while LUC 6, 8, and 10 use gamma = 0.56.
%
%   WindSpeed is an optional input to the function, but is required for the
%     Variable Deposition Model. WindSpeed is a scalar or vector value with
%     the same number of elements as Time, and must be in meters per
%     second. If WindSpeed is omitted, the value of 2 m/s is used as
%     default.
%
%   Temperature is an optional input to the function, but is required for
%     the Variable Deposition Model. Temperature is a scalar or vector
%     value with the same number of Elements as Time and must be in degrees
%     C. If Temperature is omitted, the value of 12 C is used as default.
% 
% 
% Output Parameters:
%   SR = The soiling ratio (SR) of a tilted PV panel, this is a number
%   between 0 and 1
%
% References
%   [1] M. Coello, L. Boyle. A Simple Model for Predicting Time Series
%   Soiling of Photovoltaic Panels. IEEE PVSC/ World Conference on
%   Photoltaic Energy Conversion (WCPEC) 2018 (based on pre-print of paper 
%   to be published in IEEE JPV with reference updated upon publishing).
%   [2] Atmospheric Chemistry and Physics: From Air Pollution to Climate
%   Change. J. Seinfeld and S. Pandis. Wiley and Sons 2001.
%
% See also PVL_MAKETIMESTRUCT 

p=inputParser;

p.addRequired('Time', @isstruct);
p.addRequired('Rain', @(x) all(isnumeric(x) & all(x>=0) & isvector(x)));
p.addRequired('RainThresh', @(x) all(isnumeric(x) & isscalar(x) & all(x>=0)));
p.addRequired('Tilt', @(x) all(isnumeric(x) & isvector(x)));
p.addRequired('PM2_5', @(x) all(isnumeric(x) & all(x>=0) & isvector(x)))
p.addRequired('PM10', @(x) all(isnumeric(x) & all(x>=0) & isvector(x)))
p.addOptional('ModelType', 2, @(x) all(isscalar(x)));
p.addOptional('RainAccPeriod', 1, @(x) all(isscalar(x) & x>=1 & isnumeric(x)));
p.addOptional('LUC', 8, @(x) all(isscalar(x)))
p.addOptional('WindSpeed', 2, @(x) all(isnumeric(x) & all(x>=0) & isvector(x)));
p.addOptional('Temperature', 12, @(x) all(isnumeric(x) & all(x>=0) & isvector(x)));
p.parse(Time, Rain, RainThresh, Tilt, PM2_5, PM10, varargin{:});

RainAccPeriod = p.Results.RainAccPeriod;
ModelType = p.Results.ModelType;
LUC = p.Results.LUC;
WindSpeed = p.Results.WindSpeed;
Temperature = p.Results.Temperature;

% Create a datenum from the time structure to handle duration easily
TimeAsDatenum = datenum(Time.year, Time.month, Time.day, Time.hour, ...
    Time.minute, Time.second);

% Following section accumulates the rain in each accumulation period
% To change this to a more frequent accumulation rate, create a
% new array based on TimeAsDatenum where the unit is the given unit of
% accumulation (e.g. for hourly accumulation you might use
% TimeAsDatenum*24)
RainAccAsDatenum = floor(TimeAsDatenum * 24 ./ RainAccPeriod);
% Determine the number of unique days, the first instance of that day, and
% create an index of which entries are part of that day
[RainAccTimes, UnqRainAccFrstVal, UnqRainAccIndx] = unique(RainAccAsDatenum);
% Accumulate rain readings according to the unique index of days
RainAtAccTimes = accumarray(UnqRainAccIndx, Rain);
% Create an accumulated rain vector, then populate the last entry of each
% date with the rain that fell on that day
AccumRain = zeros(size(Rain));
AccumRain(UnqRainAccFrstVal(2:end)-1) = RainAtAccTimes(1:end-1);
% Rain accumulated on the last day is placed at the last time (does not
% affect soiling calculation)
AccumRain(end) = RainAtAccTimes(end);

switch ModelType
    case 1 % Variable Deposition Velocity
        vd = pvl_depo_veloc(Temperature, WindSpeed, LUC);
    case 2 % Fixed Settling Velocity in m/s
        vd(1,2) = 0.004;
        vd(1,1) = 0.0009;
    case 3 % Fixed Deposition Velcoity in m/s
        vd(1,2) = 0.0917;
        vd(1,1) = 0.0015;
    otherwise
        warning('Unknown Model Type specified, using Fixed Setting Velocity model')
        vd(1,2) = 0.004;
        vd(1,1) = 0.0009;
end


% Corrects PM10 measurement since PM2.5 is included in measurement
PMConcentration=nan(numel(TimeAsDatenum),2); % pre-allocate with NAN

PMConcentration(:,1) = PM2_5; % fill PM2.5 data in column 1
% Since PM2.5 particulate is counted in PM10 measurements, and we need only
% the particular from 2.5 micron to 10 micron, subtract the 2.5 micron
% particulate from the 10 micron particulate
PMConcentration(:,2) = PM10 - PM2_5; % fill in PM2.5-PM10 data in column 2

% For any instances where 2.5 micron particulate was more thatn 10 micron,
% set the amount of 2.5 to 10 micron particulate to 0
PMConcentration(PM10 - PM2_5 < 0 , 2) = 0; 

% EPA concentrations reported in micrograms/m^3, we need it in g/m^3
PMConcentration = PMConcentration .* 1e-6;

% Calculate the mass deposition rate of each particulate size
F=PMConcentration .* vd; % g * m^-2 * s^-1, by particulate size
HorizontalTotalMassRate = F(:,1) + F(:,2); % g * m^-2 * s^-1, total

% Mass depositition rate on a tilted surface. In the event of a changing
% tilt (e.g. tracking system) the mass rate is adjusted by the tilt angle
% at each time step, and thus a linear interpolation of changing tilt
% angles are used (when using trapezoidal integration).
TiltedMassRate = HorizontalTotalMassRate .* cosd(Tilt);

% Mass on the tilted surface if no rain occurs. Uses trapezoidal
% integration and determines the amount of mass deposited on a tilted
% module.
TiltedMassNoRain = cumtrapz(TimeAsDatenum*86400, TiltedMassRate);

% Create a Tilted soiling mass that will account for rain cleaning
TiltedMass = TiltedMassNoRain;
% For each rain accumulation period, check to see if enough rain fell to meet or
% exceed the cleaning threshold. If the threshold was met or exceeded,
% then subtract from each subsuequent time step the mass of future time
% steps. This essentially "resets" the cumulative soiling accumulation when
% the rain meets or exceeds the stated thresshold. Note that this means
% that the mass will never be reported as zero because the cleaning is
% applied immediately after a time step, and new mass accumulates between
% the cleaning and the first time reported after the cleaning.
% Note that this could be done more clearly by going through each time step
% and accumulating or clearing based on daily rainfall, but this is much
% faster, and can utilize trapezoidal soiling integration.
for cntr1 = 1:numel(RainAtAccTimes)-1
    if RainAtAccTimes(cntr1) >= RainThresh
        TiltedMass(UnqRainAccFrstVal(cntr1+1):end) = TiltedMass(UnqRainAccFrstVal(cntr1+1):end)-TiltedMass(UnqRainAccFrstVal(cntr1+1)-1);
    end
end

% Soiling rate is the % transmittance reduction using Heagzy equation
SoilingRate = 34.37 * erf(0.17 .* (TiltedMass.^0.8473));

% Soiling Ratio is the % soil transmittance (0.95 = 5% soiling loss)
SR = (100 - SoilingRate) ./ 100;

end

% This function creates deposition velocities for the variable deposition
% model.
function vd = pvl_depo_veloc(T, WindSpeed, LUC)

    % convert temperature into Kelvin 
    T = T + 273.15;

    % save wind data
    u = WindSpeed;

    g=9.81;         %gravity in m/s^2
    Na=6.022*10^23; %avagadros number
    R=8.314;        %Universal gas consant in m3Pa/Kmol
    k=1.38*10^-23;  %Boltzmann's constant in m^2kg/sK
    P=101300;       %pressure in Pa
    rhoair= 1.2041; %density of air in kg/m3
    z0=1;
    rhop=1500;      %Assume density of particle in kg/m^3
    
    switch LUC
        case {1, 4}
            gamma = 0.56;
        case {6, 8, 10}
            gamma = 0.54;
        otherwise
            warning('Unknown Land Use Category (LUC), assuming LUC 8.');
            gamma = 0.54;
    end


    %Dimeter of particle in um
    Dpum=[2.5,10];
    Dpm=Dpum*10^-6;   %Diameter of particle in m

    %Calculations
    mu=1.8*10^-5.*(T./298).^0.85;      %viscocity of air in kg/m s
    nu=mu/rhoair;
    lambda1=2*mu./(P.*(8.*0.0288./(pi.*R.*T)).^(1/2));   %mean free path
    
    Cc = 1+2.*[lambda1 lambda1]./Dpm.*(1.257+0.4.*exp(-1.1.*Dpm./([lambda1 lambda1].*2))); %slip correction coefficient

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Calculate vs
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    vs = rhop.*Dpm.^2 .* (g .* Cc ./(mu.*18)); %particle settling velocity

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Calculate rb
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ustar = NaN(size(u)); % pre-allocate ustar
    %Equation 11.66 in Ramaswami (and 16.67 and Sienfeld &Pandis)
    ustar(u > 0) = 0.4 * u(u > 0) ./ log(10/z0);
    ustar(u<=0) = 0.001;

    D=k*T.*(Cc./(3*pi*mu*Dpm));

    Sc=nu./D;
    %gamma=.56;      %for urban
    %alpha=1.5;     %for urban      
    EB=Sc.^(-1 * gamma);
    St = vs .* (ustar .^ 2) ./ (g .* nu);
    
    EIM=10.^(-3./St);   %For smooth surfaces
    %EIM=((St)./(0.82+St)).^2;

    R1=exp(-St.^(1/2));  %percentage of particles that stick

    rb=1./(3*(EB+EIM).*ustar.*R1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Calculate ra
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    a=[-0.096,-0.037,-0.002,0,0.004,0.035];
    b=[0.029,0.029,0.018,0,-0.018,-0.036];
    
    % For wind speeds <= 3, use a = -0.037 and b = 0.029
    % For wind speeds >3 and <=5, use a = -.002, b = 0.018
    % For wind speeds > 5, use a = 0, b = 0
    avals = a(2) .* ones(size(u));
    avals(u>3) = a(3);
    avals(u>5) = a(4);
    
    bvals = b(2) .* ones(size(u));
    bvals(u>3) = b(3);
    bvals(u>5) = b(4);
    
    L = 1 ./ (avals + bvals .* log(z0));
    
    zeta0=z0./L;
    zeta=10./L;
    eta = ((1-15.*zeta).^(1./4));
    eta0 = ((1-15.*zeta0).^(1./4));
    
    ra = NaN(size(zeta)); % Preallocate memory
    ra(zeta == 0) = (1 ./ (0.4 .* ustar(zeta == 0))) .* log(10 ./ z0);
    ra(zeta > 0) = (1./(0.4.*ustar(zeta > 0))).*(log(10./z0) + ...
        4.7.*(zeta(zeta > 0)-zeta0(zeta > 0)));
    ra(zeta < 0) = (1 ./ (0.4 .* ustar(zeta < 0))) .* (log(10 ./ z0) + ...
        log((eta0(zeta < 0).^2 + 1) .* (eta0(zeta < 0)+1).^2 ./ ...
        ((eta(zeta < 0).^2 + 1) .* (eta(zeta < 0)+1).^2)) + ...
        2 .* (atan(eta(zeta < 0))-atan(eta0(zeta < 0))));


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Calculate vd and mass flux
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    vd=1./(ra+rb)+vs;

end
