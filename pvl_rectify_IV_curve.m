function [I V] = pvl_rectify_IV_curve(tI, tV, Voc, Isc)
% PVL_RECTIFY_IV_CURVE ensures that Isc and Voc are included in a IV curve
% and removes duplicate voltage and current points
%
% Syntax
%   [I V] = pvl_rectify_IV_curve(tI, tV, Voc, Isc)
%
% Description
%   pvl_rectify_IV_curve ensures that the IV curve data  
%   - increase in voltage
%   - contain no negative current or voltage values
%   - have the first data point as (0,Isc)
%   - have the last data point as (Voc,0)
%   - contain no duplicate voltage values. Where voltage values are
%   repeated, a single data point is substituted with current equal to the
%   average of current at each repeated voltage.
%
% Input:
%   tI - a vector of length N containing the current data
%   tV - a vector of length N containing the voltage data
%   Voc - a scalar containing the open circuit voltage
%   Isc - a scalar containing the short circut current
%
% Output:
%   I, V - vectors of equal length containing current and voltage,
%   respectively

% Filter out negative voltage and current values
u=tV<Voc & tV>0 & tI>0;

% Add in Voc and Isc
V=[0; tV(u); Voc];
I=[Isc; tI(u); 0];

% remove duplicate Voltage points
[C,ia,ic] = unique(V);
c = hist(ic,1:length(ia));
cfil = c>1;
ind = find(cfil);
if ~isempty(ind)
    for k=1:length(ind)
       nV = V(ia);
       nI = I(ia);
       dupfil = ic==ind(k); %identify duplicates
       nI(ind(k)) = mean(I(dupfil)); %assign mean current
    end
    V = nV;
    I = nI;
end;

end
