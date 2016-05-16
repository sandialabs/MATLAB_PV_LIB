%% pvl_rectify_IV_curve
% Ensures that Isc and Voc are included in a IV curve and removes duplicate
% voltage and current points.
%
%% Syntax
% |[I V] = pvl_rectify_IV_curve(tI, tV, Voc, Isc)|
%
%% Description
% |pvl_rectify_IV_curve| ensures that the IV curve data:
%
% * increase in voltage,
% * contain no negative current or voltage values,
% * have the first data point as (0,Isc),
% * have the last data point as (Voc,0),
% * contain no duplicate voltage values.
%
% Where voltage values are
% repeated, a single data point is substituted with current equal to the
% average of current at each repeated voltage.
%
%% Inputs
% * *|tI|* - a vector of length N containing the current data.
% * *|tV|* - a vector of length N containing the voltage data.
% * *|Voc|* - a scalar containing the open circuit voltage.
% * *|Isc|* - a scalar containing the short circut current.
%
%% Output
% * *|I, V|* - vectors of equal length containing current and voltage,
%   respectively.
%

%% Copyright 2014 Sandia National Laboratories

