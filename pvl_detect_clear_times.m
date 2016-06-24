function [clearSamples, csGHI, alpha] = pvl_detect_clear_times(GHI, Time, UTCoffset, Location, win_length, sample_interval)
% PVL_DETECT_CLEAR_TIMES identifies times with irradiance consistent with 
% clear sky conditions in a time series of GHI.
%
% Syntax
%   [clearSamples, csGHI, alpha] = pvl_detect_clear_times(GHI, Time, UTCoffset, Location, win_length, sample_interval)
%
% Description
%   Detects clear sky times by comparing statistics for a regular GHI 
%   time series to the Ineichen clear sky model.  Statistics are calculated
%   using a sliding time window (e.g., 10 minutes). An iterative algorithm
%   identifies clear periods, uses the identified periods to estimate bias
%   in the clear sky model and/or GHI data, adjusts the clear sky model and
%   repeats.  Code handles GHI data with some irregularities, i.e., missing
%   values or unequal data spacing.  Code execution can be made
%   significantly faster if equally spaced and complete data can be
%   assumed.
%   Clear times are identified by meeting 5 criteria, thresholds
%   for which are hardcoded in this version.  Values for these thresholds
%   are appropriate for 10 minute windows of 1 minute GHI data.
%
% Inputs:   
%   GHI - a Nx1 vector of GHI values (W/m2)
%   Time - a Nx1 vector of datenum with equal time spacing
%   UTCoffset - scalar UTC offset (e.g. EST = -5)
%   Location - a struct with the following elements, note that all
%      elements must be scalars in this application.
%   Location.latitude - scalar latitude in decimal degrees (positive is
%      northern hemisphere)
%   Location.longitude - scalar longitude in decimal degrees (positive is
%      east of prime meridian)
%   Location.altitude - scalar height above sea level in meters.
%      While altitude is optional in many uses, it is required in this
%      model implementation.
%   win_length - length of sliding time window in minutes, integer
%   sample_interval - nominal minutes between each GHI sample time, integer
%
% Output:   
%   clearSamples - column vector of logical values with True indicating 
%    a clear sample in the input GHI time series
%   csGHI - column vector with scaled clear sky GHI derived from the
%   Ineichen model
%   alpha - scaling factor applied to Ineichen model to obtain output csGHI
%   
% References
%   [1] Reno, M.J. and C.W. Hansen, "Identification of periods of clear sky
%   irradiance in time series of GHI measurements" Renewable Energy, v90, 
%   p. 520-531, 2016.
%
% Notes:
%   Initial implementation by Matthew Reno. Modifications for computational
%   efficiency by Joshua Patrick and Curtis Martin.
%


% Threshold values for each criterion.  Values are appropriate for 10
% minute windows with 1 minute intervals between GHI data.

% MeanDiff: threshold value in W/m2 for first criterion, i.e., agreement between mean 
% values of GHI in each interval, see Eq. 6 in [1]
MeanDiff = 75;  

% MaxDiff: threshold value in W/m2 for second criterion, i.e., agreement
% between maxima of GHI values in each interval, see Eq. 7 in [1]
MaxDiff = 75;

% LineLengthLower and LineLengthUpper: threshold values (unit is weird) for
% third criterion, i.e., agreement between line lengths of GHI data and 
% clear sky GHI, see Eq. 8 in [1].  Criterion is satisfied when line length
% LL meets LineLengthLower <= LL <= LineLengthUpper
LineLengthLower = -5; 
LineLengthUpper = 10;

% VarDiff: threshold value in 1/seconds for the fourth criterion, i.e.,
% agreement between normalized standard deviations of rate of change in
% irradiance, see Eq. 9 through Eq. 11 in [1]
VarDiff = 0.005;

% SlopeDev: threshold value in W/m2 for the fifth criterion, i.e.,
% agreement between largest magnitude of change in successive GHI values,
% see Eq. 12 through 14 in [1]
SlopeDev = 8;

dv = datevec(Time);
dates = datenum(dv(:,1:3));
days = unique(dates);
ivec = (1:length(Time))';

% Get standard Ineichen clear sky GHI
csGHI0 = pvl_clearsky_ineichen(pvl_maketimestruct(Time,UTCoffset),Location);
% replace NaNs with 0
csNaNs = isnan(csGHI0);
csGHI0(csNaNs) = 0;

% Create index matrix for windows

max_samples_per_window = 2 * win_length / sample_interval;  % upper bound on number of data samples in an interval
min_samples_per_window = floor(0.8 * win_length / sample_interval) ; % lower bound on number of data samples in an interval

% initialize some arrays
clearMean = NaN(size(GHI));
clearMax = clearMean;
clearMaxSlope = clearMean;
timeDiff = zeros(max_samples_per_window, length(clearMean));
clearGHIDiff = timeDiff;

measuredMean  = clearMean;
measuredMax  = clearMean;
measuredMaxSlope  = clearMean;
measuredSqSlope  = clearMean;
measuredStdSlope  = clearMean;
measuredLineLength = clearMean;


k = 0; % counts number of intervals meeting data size criteria, i.e., samples between min_samples_per_window and max_samples_per_window
for j = 1:length(Time)
    tu = Time>=Time(j) & Time<(Time(j)+ (win_length + 1)/1440);
    if sum(tu) >= min_samples_per_window && sum(tu) <= max_samples_per_window  % at least 8 measurements within a window but not more than 20
        k = k + 1;
        tidx = ivec(tu);
        w = diff(Time(tidx)*1440);  % time intervals in minutes between samples
        
        % Calculate parameters (all but line length) for alpha=1 case
        wx = [w; 0]/sum(w); % leave last point out of mean
        clearMean(j) = sum(wx.*csGHI0(tidx));
        clearMax(j) = max(csGHI0(tidx));
        clearMaxSlope(j) = max(abs(diff(csGHI0(tidx),1,1)./w));
        tmp_diffGHI = diff(csGHI0(tidx));
        clearGHIDiff(1:length(tmp_diffGHI),j) = tmp_diffGHI;
        timeDiff(1:length(tmp_diffGHI),j) = w;

        % Calculate parameters for measured GHI
        measuredMean(j) = sum(wx.*GHI(tidx));
        measuredMax(j) = max(GHI(tidx));
        tmpSlope = diff(GHI(tidx),1,1)./w;
        measuredMaxSlope(j) = max(abs(tmpSlope));
        measuredSqSlope(j) = sum(tmpSlope.^2);
        measuredStdSlope(j) = std(tmpSlope);
        measuredLineLength(j) = sum(sqrt(diff(GHI(tidx),1,1).^2 + w.^2)); 
    else
        % too much missing data or too many data samples
        disp(['Data problem between ' datestr(min(Time(tu))) ' and ' datestr(max(Time(tu))) ...
            ' have ' num2str(sum(tu)) ' values']);
    end
end

% Normalized variance of measured GHI < 0.005 criterion doesn't depend on
% clear sky GHI
c4 = measuredStdSlope ./ measuredMean < VarDiff;

% initialize state for iterative process
alpha = 1;      %start off with Ineichen clear sky model 

% Start Iterative process of finding clear days, fitting clear sky model, finding clear days, ....
% Continue iterations until not finding any more clear days or 20 iterations
for ii=1:20

    % Update clear sky model by scaling Ineichen
    csGHI = alpha * csGHI0;
    
    % Adjust Clear Sky Model parameters for scaling by alpha
    clearLineLength = sum(sqrt((alpha*clearGHIDiff).^2 + timeDiff.^2 ));
    clearLineLength = clearLineLength(:);

    % Evaluate comparison criteria
    c1 = abs(measuredMean - alpha * clearMean) < MeanDiff;
    c2 = abs(measuredMax - alpha * clearMax) < MaxDiff;
    c3 = (LineLengthLower < measuredLineLength - clearLineLength) & (measuredLineLength - clearLineLength < LineLengthUpper) ;
    c5 = (measuredMaxSlope - alpha * clearMaxSlope) < SlopeDev;
    c6 = clearMean ~= 0 & ~isnan(clearMean);            % window includes some daylight
    
    % Identify windows meeting all five (six) criteria
    clearWindows = c1 & c2 & c3 & c4 & c5 & c6;

    % Daily clearness is proportion of samples that are clear
    dailyClearness = histcounts(Time(clearWindows),[days]) ./ ...
        histcounts(Time(csGHI>25),[days]);
    
    % Save current state
    previousAlpha = alpha;

    % Get new alpha by adjusting clear sky model to match clear times
    % Need to minimize RMSE between GHI and alpha * csGHI0
    
    clearSamples = false(size(GHI));

    for j = 1:length(Time)
        if clearWindows(j)
            tu = Time>=Time(j) & Time<(Time(j)+win_length/1440);
            clearSamples(tu) = true;
        end
    end

    RMSE = @(x)sqrt(mean((GHI(clearSamples) - x*csGHI0(clearSamples)).^2));
    alpha = fminsearch(RMSE, alpha);  % update scaling factor on clear sky model

    if round(alpha*10000)==round(previousAlpha*10000)
        previousAlpha = alpha; %#ok<NASGU>
        break
    end
        
end  % for loop

csGHI = alpha * csGHI0;

end  % function


