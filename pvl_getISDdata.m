function outfilen = pvl_getISDdata(lat,lon,yr,archive)
%
% getISDdata: fetch data for one year of integrated surface data for location
%
% Syntax
%   data = getISDdata(latitude,longitude,year,archive)
%
% Inputs:
%   latitude  - latitude of location of interest
%   longitude - longitude of location of interest
%   year      - desired year
%   archive   - string specifying path to local ISD archive.  Current 
%               directory is used if archive is not specified.
%
% Output:
%   outfilen - name for the returned data file
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
    display('Function pvl_getISDdata requires R2013b or later')
    return
end
    
%% Get index data
% Change pwd to a local directory if you want to use an archive
if nargin < 4
    archive = pwd;
end

isdfile = 'isd-history.csv';
fname = [archive '\' isdfile];
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
index.LAT(strcmp('',index.LAT)) = {'99999'};
index.LON(strcmp('',index.LON)) = {'99999'};
index.ELEV(strcmp('',index.ELEV)) = {'99999'};
index.LAT = str2num(char(index.LAT)); %#ok<*ST2NM>
index.LON = str2num(char(index.LON));
index.ELEV = str2num(char(index.ELEV));
index.BEGIN = datenum(char(index.BEGIN),'yyyymmdd');
index.END = datenum(char(index.END),'yyyymmdd');

%% Find row for closest location (might iterate through several locations)
dist = hypot(index.LAT-lat,index.LON-lon);
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
fname = [archive '\' outfilen];
if ~exist(fname,'file')
    addr = sprintf('ftp://ftp.ncdc.noaa.gov/pub/data/noaa/%4d/',yr);
    gunzip([addr '/' outfilen '.gz'],archive);
    %% Alt approach if gunzip doesn't seem to be working
    % % cd(ftphandle,num2str(yr));
    % % gunzip(mget(ftphandle,[fname '.gz']));
    % % dat = readISH(fname);
    % % close(ftphandle);
end

