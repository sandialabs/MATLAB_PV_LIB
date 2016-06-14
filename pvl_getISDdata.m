function outfilen = pvl_getISDdata(latitude,longitude,year,archive)
%
% getISDdata: fetch data for one year of integrated surface data for location
%
% Syntax
%   data = getISDdata(latitude,longitude,year,archive)
%
% Inputs:
%   latitude  - latitude of location of interest in decimal degrees
%   longitude - longitude of location of interest in decimal degrees
%   year      - desired year
%   archive   - string specifying path to local ISD archive.  Current 
%               directory is used if archive is not specified.
%
% Output:
%   outfilen - name for the returned data file in the form "USAF-WBAN-yyyy"
%      where USAF is the USAF ID of the station and WBAN is the WBAN ID of 
%      the system, as provided by the ISD.
%
% Notes:
%   This function copies up to two files from NOAA's ftp site to the
%   current directory.  One, isd-history.csv, is a listing of all sites in
%   the database.  The other is the file containing the data requested.
%   These two files can be re-directed to an archive folder by supplying
%   the path to the archive as input archive.  If the files already exist
%   in the location pointed to by archive, the local files are used.
%
%   Code requires Matlab 2013b or later.
%
% References:
%   National Climatic Data Center, "Federal Climate Complex Data
%   Documentation for Integrated Surface Data 21 May 2014," 
%   ftp://ftp.ncdc.noaa.gov/pub/data/noaa/ish-format-document.pdf,
%   accessed 5 June 2014.

%% check Matlab version

vers = version;

pos = regexpi(vers,'(R');
versyr = str2num(vers(pos+2:pos+5));
mod = vers(pos+6);
if versyr<2013 || ((versyr==2013) && mod=='a')
    display(['Matlab ' vers ' detected'])
    error('Function pvl_getISDdata requires R2013b or later')
end
%% Check the input data
p = inputParser;
p.addRequired('latitude',@(x) all(isscalar(x) & isnumeric(x) & x<=90 & x>=-90));
p.addRequired('longitude', @(x) all(isscalar(x) & isnumeric(x) & x<=180 & x>=-180));
p.addRequired('year',@(x) all(isscalar(x) & isnumeric(x) & (x >= 1901)));
p.addOptional('archive', @(x) ischar(x));
p.parse(latitude, longitude, year, varargin{:});

defaultchecker = {'archive'};
if ~any(strcmp(defaultchecker,p.UsingDefaults))
    archive = pwd;
end

lat = p.Results.latitude;
lon = p.Results.longitude;
yr = p.Results.year;

%% Get index data
% Change pwd to a local directory if you want to use an archive

isdfile = 'isd-history.csv';
fname = [archive filesep isdfile];
if ~exist(fname,'file')
    % Get index file
    ftphandle = ftp('ftp.ncdc.noaa.gov');
    cd(ftphandle,'pub/data/noaa');
    mget(ftphandle,isdfile,archive);
    close(ftphandle);       % will reconnect automatically if we need this later
end
%%
% Read/parse index
warning('off'); % Eliminates warning that variable names changed
index = readtable(char(fname));
warning('on');

%%
index.Properties.VariableNames{'ELEV_M_'} = 'ELEV';
index.USAF = char(index.USAF);
index.WBAN = char(index.WBAN);
index.LAT(strcmp('',index.LAT)) = {'NaN'};
index.LON(strcmp('',index.LON)) = {'NaN'};
index.ELEV(strcmp('',index.ELEV)) = {'NaN'};
index.LAT = str2num(char(index.LAT)); %#ok<*ST2NM>
index.LON = str2num(char(index.LON));
index.ELEV = str2num(char(index.ELEV));
index.BEGIN = datenum(char(index.BEGIN),'yyyymmdd');
index.END = datenum(char(index.END),'yyyymmdd');

%% Find row for closest location (might iterate through several locations) using haversine formula
% calculate the argument for the haversine formula, make sure the argument
% is between -1 and 1
arg = sqrt((sind((index.LAT - lat) ./ 2)).^2 + cosd(index.LAT).*cosd(lat).*(sind((index.LON-lon)./2)).^2);
arg(arg < -1) = -1;
arg(arg > 1) = 1;
% compute the distance in meters between the provided lat/lon and the each lat/lon from the isd.
dist = 2 .* 6371008.8 .* asin(arg); % mean earth radius is 6371.0088 km
[~,idx] = sort(dist);

% Find first row that contains year (might iterate through several years)
ii = find(datenum(yr,1,1)>=index.BEGIN(idx) & datenum(yr,12,31)<=index.END(idx),1);

if isempty(ii)
    % I don't think this will happen.  The nearest station may be far away,
    % but should be able to find one.
    error('No station found')
end

row = idx(ii);
%%

% form file name
outfilen = sprintf('%s-%s-%4d',index.USAF(row,:),index.WBAN(row,:),yr);

% 
% get file from ftp server if not in archive
fname = [archive filesep outfilen];
if ~exist(fname,'file')
    addr = sprintf('ftp://ftp.ncdc.noaa.gov/pub/data/noaa/%4d/',yr);
    try
        gunzip([addr '/' outfilen '.gz'],archive);
    catch
        error(['Unable to reach the ftp site for the integrated surface database.' ...
            ' The data file cannot be accessed.']);
    end
    %% Alt approach if gunzip doesn't seem to be working
    % % cd(ftphandle,num2str(yr));
    % % gunzip(mget(ftphandle,[fname '.gz']));
    % % dat = readISH(fname);
    % % close(ftphandle);
end

