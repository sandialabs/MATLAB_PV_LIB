%% pvl_getISDdata
% Fetch one year of data from the Integrated Surface Database.
%
%% Syntax
% |fname = getISDdata(latitude,longitude,year,archive)|
%
%% Inputs
% * *|latitude|*  - latitude of location of interest.
% * *|longitude|* - longitude of location of interest (negative for west of prime meridian).
% * *|year|*      - desired year.
% * *|archive|*   - string specifying path to local ISD archive.  The current 
%               directory is used if archive is not specified.
%
%% Output
% * *|fname|* - name of the file containing returned data.
%
%% Description
% This function copies two files from NOAA's ftp site <'ftp://ftp.ncdc.noaa.gov/pub/data/noaa/> to the
% current directory.  The file |isd-history.csv| is a listing of all sites in
% the database.  The second file contains one year of data from the station closest to the requested location.
%
% These two files can be re-directed to an archive folder by supplying
% the path to the archive as |archive|.  If the files already exist
% in the location pointed to by |archive|, the local files are used. 
% Use <pvl_readISH_help.html |pvl_readISH|> to read the ISD data files.
%
% Requires Matlab 2013b or later.
%
%% Example
archive = '..\Example Data';
fname = pvl_getISDdata(35.11,-106.61,2010,archive); % get weather data for Albuquerque, NM
dir([archive '\isd-history.csv'])  % index of all data in ISD
dir([archive '\' fname]) % data for 2010 for Albuquerque, NM
data = pvl_readISH([archive '\' fname]);
data(1,1:18)

%% References
% [1] National Climatic Data Center, Federal Climate Complex Data
% Documentation for Integrated Surface Data 21 May 2014,
% ftp://ftp.ncdc.noaa.gov/pub/data/noaa/ish-format-document.pdf,
% accessed 5 June 2014. 

%% See Also 
% <pvl_readISH_help.html |pvl_readISH|> 
%%
% Copyright 2015 Sandia National Laboratories