%% pvl_readISH
% Read text file of Integrated Surface History (ISH) data.
%
%% Syntax
% |data = readISH(ISHname)|
%
%% Description
% Reads an ASCII file containing one year of from the Integrated Surface
% Database <ftp://ftp.ncdc.noaa.gov/pub/data/noaa/> for a single station. Data are copied
% into a table with columns for each data field.
%
% Requires Matlab 2013b or later.
%
%% Input
% * *|ISHname|* - string specifying the file (including path) containing the ISH data.
%
%% Output
% * *|data|* - table containing ISH data, one row per observation,
%         one column for each mandatory data element. The
%         last column in the table contains a string with the
%         additional data and remarks.
%
% Data columns are as follows (see [1] for definition of each field):
%
% <html>
% <table border=1>
% <tr><td>Label</td>      <td>Type</td>      <td>Description (units)</td></tr>
% <tr><td>len</td>        <td>int</td>       <td>length of variable portion of the record</td></tr>
% <tr><td>catno</td>      <td>string</td>    <td>AF Catalog Station Number</td></tr>
% <tr><td>wban</td>       <td>string</td>    <td>WBAN number</td></tr>
% <tr><td>date</td>       <td>string</td>    <td>Date of observation, in YYYYMMDD format</td></tr>
% <tr><td>time</td>       <td>string</td>    <td>Time of observation, in HHMM format</td></tr>
% <tr><td>source</td>     <td>char</td>      <td>Data source code</td></tr>
% <tr><td>lat</td>        <td>double</td>    <td>Observation latitude (degrees North)</td></tr>
% <tr><td>long</td>       <td>double</td>    <td>Observation longitude (degrees East)</td></tr>
% <tr><td>type</td>       <td>string</td>    <td>Report type code</td></tr>
% <tr><td>elev</td>       <td>double</td>    <td>Observation elevation (meters)</td></tr>
% <tr><td>station</td>    <td>string</td>    <td>Station call letters</td></tr>
% <tr><td>qc</td>         <td>string</td>    <td>Name of quality control process applied</td></tr>
% <tr><td>winddir</td>    <td>double</td>    <td>Wind direction (degrees clockwise from N)</td></tr>
% <tr><td>winddirq</td>   <td>char</td>      <td>Wind direction quality code</td></tr>
% <tr><td>windtype</td>   <td>char</td>      <td>Wind observation type code</td></tr>
% <tr><td>windspeed</td>  <td>double</td>    <td>Wind speed (meters/second)</td></tr>
% <tr><td>windspeedq</td> <td>char</td>      <td>Wind speed quality code</td></tr>
% <tr><td>ceil</td>       <td>double</td>    <td>Cloud ceiling height (meters AGL)</td></tr>
% <tr><td>ceilq</td>      <td>char</td>      <td>Cloud ceiling quality code</td></tr>
% <tr><td>ceilmethod</td> <td>char</td>      <td>Cloud ceiling measurement method code</td></tr>
% <tr><td>CAVOK</td>      <td>char</td>      <td>Ceiling and Visibility OK code</td></tr>
% <tr><td>visdist</td>    <td>double</td>    <td>Visibility distance (meters)</td></tr>
% <tr><td>visdistq</td>   <td>char</td>      <td>Visibility distance quality code</td></tr>
% <tr><td>visvar</td>     <td>char</td>      <td>Visibility variability code</td></tr>
% <tr><td>visvarq</td>    <td>char</td>      <td>Visibility variability quality code</td></tr>
% <tr><td>temp</td>       <td>double</td>    <td>Air temperature (degrees C)</td></tr>
% <tr><td>tempq</td>      <td>char</td>      <td>Air temperature quality code</td></tr>
% <tr><td>dewpt</td>      <td>double</td>    <td>Dew point temperature (degrees C)</td></tr>
% <tr><td>dewptq</td>     <td>char</td>      <td>Dew point temperature quality code</td></tr>
% <tr><td>pressure</td>   <td>double</td>    <td>Barometric pressure (hectopascals)</td></tr>
% <tr><td>pressureq</td>  <td>char</td>      <td>Barometric pressure quality</td></tr>
% <tr><td>remainder</td>  <td>string</td>    <td>Optional part of record (ADD + REM)</td></tr></table>
% </html>
%

%% Example
archive = '..\Example Data';
fname = pvl_getISDdata(35.11,-106.61,2010,archive); % get weather data for Albuquerque, NM
dir([archive '\isd-history.csv'])  % index of all data in ISD
dir([archive '\' fname]) % data for 2010 for Albuquerque, NM
data = pvl_readISH([archive '\' fname]);
data(1,1:18)

%% Reference
% [1] National Climatic Data Center, Federal Climate Complex Data
% Documentation for Integrated Surface Data 21 May 2014,
% ftp://ftp.ncdc.noaa.gov/pub/data/noaa/ish-format-document.pdf,
% accessed 5 June 2014. 
%% See Also 
% <pvl_getISDdata_help.html |pvl_getISDdata|> 
%%
% Copyright 2015 Sandia National Laboratories