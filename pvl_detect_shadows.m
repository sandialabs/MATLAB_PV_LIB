function out = pvl_detect_shadows(Time, GHI, interval, site_info, dbg)
% pvl_detect_shadows identifies shading on a GHI instrument from local fixed
% structures such as wires and poles.
%
% Syntax
%   out = pvl_detect_shadows(Time, GHI, dbg)
%
% Description
%   pvl_detect_shadows uses morphological image processing methods to
%   identify shadows from fixed local objects in GHI data. GHI data are 
%   are assumed to be reasonably complete with relatively few missing values
%   and at a fixed time interval nominally of 1 minute over the course of 
%   several months. Shadows are assumed to be relatively short duration. 
%   The algorithm forms a 2D image of the GHI data by arranging time of day 
%   along the x-axis and day of year along the y-axis. Rapid change in GHI
%   in the x-direction is used to identify edges of shadows; continuity in
%   the y-direction is used to separate local object shading from cloud
%   shadows.
%
% Inputs
%   Time - a Kx1 vector of datetimes in UTC corresponding to the measured GHI
%   GHI - a Kx1 vector of measured GHI (W/m2) 
%   interval - an integer number of minutes between data points
%   site_info - a structure with the following fields:
%     lat - site latitude
%     long - site longitude, positive is east of prime meridian
%     alt - site altitude in meters
%     UTCoffset - in hours, positive for east of prime meridian
%   dbg - a boolean; additional plots are created if set to True
%
% Output
%   out - a structure with the following fields:
%     localtime - a Kx1 vector of datetime in adjusted by UTCoffset
%     GHI - copy of input GHI in M x N = (samples per day) x (days) format
%     boosted - M x N array of GHI after adjustment for dynamic range
%     alpha - scalar used to adjust dynamic range
%     gradient - M x N array of morphological gradient of boosted
%     clouds - M x N array of boolean with True indicating non-clearsky conditions
%     night - M x N array of boolean with True indicating clearsky GHI<100 W/m2
%     thresholded - M x N array of boolean with True indicating a large morphological gradient 
%     wires - M x N array of boolean with True indicating a shadow
% 
% Notes
%   Morphological image functions are defined as methods in class GHIImage.
%   Method names are common in literature, and the Matlab Image Processing
%   toolbox uses the common naming convention.  To avoid namespace conflict
%   with the Matlab toolbox, the PVLib implementation uses a class
%   definition.
%
% References
%   [1] Martin, C. E., Hansen, C. W., An Image Processing Algorithm to 
%   Identify Near-Field Shading in Irradiance Measurements, preprint 2016
%   [2] Reno, M.J. and C.W. Hansen, "Identification of periods of clear sky
%   irradiance in time series of GHI measurements" Renewable Energy, v90, 
%   p. 520-531, 2016.

if ~mod(1440, interval)==0
    fprintf('Specified data interval is not an integer evenly dividing 1440 (minutes in a day)')
    return
end

if nargin < 5
    dbg = false;
end


M = 1440/interval;   % Number of data points in a day = rows in the image
N = floor(length(GHI)/M);  % number of whole days
% Drop partial day at the end, if present
Time(M*N+1:end) = [];
GHI(M*N+1:end) = [];

% Get clear sky GHI
fprintf('\nComputing clear sky estimate...');
Location = pvl_makelocationstruct(site_info.lat,site_info.long,site_info.alt);
[csGHI,~,~] = pvl_clearsky_ineichen(pvl_maketimestruct(Time,0),Location);

% Fill in any NaNs in GHI by linear interpolation in time
idx = isnan(GHI);
GHI(idx) = interp1(find(~idx),GHI(~idx),find(idx));

% Reshape data into image format

tGHI = GHIImage(reshape(GHI,M,N));  % create instance to get access to methods for image processing

% Scale smoothed GHI, csGHI to have peak of 1000. Removes any small bias in
% the clear sky model.
peakGHI = max(max(smooth2a(tGHI,15,15)));

sGHI = 1000 * GHI / peakGHI;
scsGHI = 1000 * csGHI / max(csGHI(:));

% Find clouds using Reno method [2] on a 50 minute window
fprintf('\nRemoving clouds...');
out.clouds = find_clouds(sGHI, Time, scsGHI, 50);

% Reshape data into image format
tGHI = GHIImage(reshape(sGHI,M,N));
scsGHI = reshape(scsGHI,M,N);
out.clouds = reshape(out.clouds,M,N);
out.night = scsGHI< 100 | isnan(scsGHI);

% capture data so far
out.localtime = Time + site_info.UTCoffset/24;
out.GHI = reshape(GHI,M,N);

% % % Interpolate across days to replace clouds
OK = ~out.night & ~out.clouds;
yy = 1:N;
for ii = 1:M
    idx = OK(ii,:);
    if sum(idx)>2
        tGHI.Value(ii,~idx) = interp1(yy(idx),tGHI.Value(ii,idx),yy(~idx));
    end
end
tGHI.Value(out.night) = nan;

% Normalize the GHI and dampen the dynamic range where the clear sky model 
% may have large errors, e.g., at very low sun elevation
alpha = 2000;
boosted = 1000 * (tGHI.Value+alpha) ./ (scsGHI+alpha);

out.boosted = boosted;
out.alpha = alpha;


% Perform morphological processing to locate wires
% THIN WIRES APPROACH
% Threshold for GHI/(csGHI+alpha): 15
% Threshold for (GHI+alpha)/(csGHI+alpha): 3
mask = ones(3,1); % use 3 samples on same day (x=3, y=1) for morphological gradient mask
mingg = 2; % threshold value for morphological gradient
A = GHIImage(boosted).graygradient(mask);
BWA = A > mingg; % filter for large gradients
connectivity = 8; % can only be 4 for vertical and horizontal neighbors, or 8 for all neighbors
tmp = GHIImage(BWA).bwclose(ones(3,1));
wiresA = GHIImage(tmp).bwareaopen(connectivity, 200); % find connected shadows comprising at least 200 pixels using all neighbors (8)

% remove clouds connected to wire shadows
wiresB = clean_wires(wiresA);

if dbg
    bw = -1; % colormap
    f = make_shadow_figure((out.GHI+1).*wiresA, out.localtime, bw);
    title('Detected shadows before cleaning');
    f = make_shadow_figure((out.GHI+1).*wiresB, out.localtime, bw);
    title('Detected shadows after cleaning');
end


%% Developmental code to try different structuring elements for the gradient 
% % % MEDIUM WIRES APPROACH
% % %$$$ Need to figure out threshold & post processing
% % % Threshold for GHI/(csGHI+alpha): 20
% % % Threshold for (GHI+alpha)/(csGHI+alpha): 5
% % mask = ones(7,1); mingg = 4;
% % B = graygradient(residual,mask);
% % BWB = B > mingg;
% % wiresB = bwareaopen(BWB,8,200);  
% % % wiresB = bwaspectopen(wiresB,8,.85);
% % 
% % % FAT WIRES APPROACH
% % %$$$ leaves some blobs in SMUD15; C4 is same as C8 for SMUD15
% % %$$$ also probably need to erode after getting rid of blobs 
% % % Threshold for GHI/(csGHI+alpha): 20
% % % Threshold for (GHI+alpha)/(csGHI+alpha): 25
% % mask = ones(15,1); mingg = 5;
% % C = graygradient(residual,mask);
% % BWC = C > mingg;
% % % wires = bwareaopen(BW,8,800);  
% % wiresC = bwareaopen(bwopen(BWC, ones(1,3)),8,800);  
% % % wiresC = bwaspectopen(wiresC,8,.8);

fprintf('\n');

out.gradient = A;
out.thresholded = BWA;
out.wires = wiresB;

%$$$@ Alternative "fat wire" masks
% "Cross" looks at same minute both forward and back one day
% mask = zeros(15,3); mask(:,2) = 1; mask(8,:) = 1;
% Look only back one day
% mask = zeros(15,3); mask(:,2) = 1; mask(8,3) = 1;
% Look only forward one day: mask(8,1) = 1;
% mask = zeros(15,3); mask(:,2) = 1; mask(8,1) = 1;
% Look at minute before and after on day before and after
% mask = zeros(15,3); mask(:,2) = 1; mask(7:9,:) = 1;
% Look at same minute back 2 days:
% mask = zeros(15,5); mask(:,3) = 1; mask(8,5) = 1;
% Look at 3 minutes back 2 days
% mask = zeros(15,5); mask(:,3) = 1; mask(7:9,5) = 1;

end


function clouds = find_clouds(GHI, Time, csGHI, window)

% private implementation of pvl_detect_clear_times.m, the algorithm to find
% cloud shadows in GHI data. Data are assumed to be complete, with no 
% missing samples, and at regular intervals. 
%
% Inputs
%   GHI - a Nx1 vector of GHI
%   Time - a Nx1 vector of datetime
%   csGHI - the clearsky GHI at each Time
%   window - the length of a window considered to distinguish clear from 
% cloudy conditions.
%
% Output
%   clouds - 
%
% References
%   [1] Reno, M.J. and C.W. Hansen, "Identification of periods of clear sky
%   irradiance in time series of GHI measurements" Renewable Energy, v90, 
%   p. 520-531, 2016.


% Need instantaneous differences (keep same lengths)
dGHI = [diff(GHI); 0];
dcsGHI = [diff(csGHI); 0];
dT = 24*60*[diff(Time); 0];    % use units of minutes
dL = hypot(dGHI,dT);
dcsL = hypot(dcsGHI,dT);

% Replace nans in csGHI with zeros
missing = isnan(csGHI);
csGHI(missing) = 0;

% Compute metrics for sliding window of 10 samples
% (window is centered on time--4 samples before, present, 5 samples after)

slidewin = ones(10,1);

% Averages
avgGHI = conv(GHI,slidewin,'same')./conv(ones(size(GHI)),slidewin,'same');
avgcsGHI = conv(csGHI,slidewin,'same')./conv(1*~missing,slidewin,'same');

% Sums
L = conv(dL,slidewin,'same');
csL = conv(dcsL,slidewin,'same');

% Variance of slope (note: assumes zero padding; inaccurate for first 4 and last 5 samples)
varwin = ones(9,1)/9;
varslope = (conv(dGHI.^2,varwin,'same') - conv(dGHI,varwin,'same').^2)./avgGHI;

% Maximum quantities
idx = min(max(hankel(1:10,10:length(GHI)+9)-4,1),length(GHI));
maxGHI = max(GHI(idx))';
maxcsGHI = max(csGHI(idx))';
maxdslope = max(abs(dGHI(idx)-dcsGHI(idx)))';

%  Apply thresholds to determine cloudy times

m1 = abs(avgGHI-avgcsGHI) < 75;
m2 = abs(maxGHI-maxcsGHI) < 75;
m3 = (-5 < (L-csL)) & ((L-csL) < 10);
% m3 = (-10 < (L-csL)) & ((L-csL) < 30);
% m4 = sigma < 0.02;
m4 = varslope < 0.005;
m5 = maxdslope < 8;

clouds = conv(1*(m1 & m2 & m3 & m4 & m5),ones(window,1),'same') == 0;

clouds = clouds & ~missing;

end


function out = clean_wires(in)

% Eliminate up to 3-day "pillars"
%  0 x x x 0
%  0 x 1 x 0
%  0 x 1 x 0
J = [0 0 0 0 0; 0 0 1 0 0; 0 0 1 0 0]; 
K = [1 0 0 0 1; 1 0 0 0 1; 1 0 0 0 1];
A1 = in & ~GHIImage(in).bwhitmiss(J,K);
A2 = in & ~GHIImage(in).bwhitmiss(flipud(J),flipud(K));
out = A1 & A2;

% Fill in any 1-minute gaps
out = out | GHIImage(out).bwhitmiss([1 0 1]');

% Eliminate 1-day spikes
out = out & ~GHIImage(out).bwhitmiss([0 1 0]);
% Fill in any 1-minute gaps
out = out | GHIImage(out).bwhitmiss([1 0 1]');

% Restore 2-minute and 1-minute gaps
out = out | GHIImage(out).bwhitmiss([0 1 0 0 1]',[0 0 1 1 0]');
out = out | GHIImage(out).bwhitmiss([1 0 1]');

% Eliminate blobs shorter than 20 days
out = GHIImage(out).bwlengthopen(4,20);

% Keep only shadows that fall within structuring elements that are
% bar-shaped and not within-day (angle = 90)
for angle = 0:5:85

    if angle < 75 
        dim = 21;
        mid = 11;
    else
        dim = 11 + 2*angle/5;
        mid = 6 + angle/5;
    end
    
%%
    if angle <= 45
%         dim = 21 + 2*angle/5;
%         mid = 11 + angle/5;
%         dim = 21;
%         mid = 11;
        cc = 1:dim;
        rr = mid + round((cc-mid)*tand(angle));
        SE = zeros(dim);
        SE(sub2ind(size(SE),rr,cc)) = 1;
    else
%         dim = 21 + 2*angle/5;
%         mid = 11 + angle/5;
        cc = 1:dim;
        rr = mid + round((cc-mid)*tand(90-angle));
        SE = zeros(dim);  
        SE(sub2ind(size(SE),rr,cc)) = 1;
        SE = SE';
    end
%%
    out = out | GHIImage(in).bwopen(SE) | GHIImage(in).bwopen(flipud(SE));
        
end

% Eliminate blobs shorter than 15 days

out = GHIImage(out).bwlengthopen(8,15);
end


