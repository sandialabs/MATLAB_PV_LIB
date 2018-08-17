function [Rs, errest] = pvl_est_Rs_Pysch(IVCurves, delI, varargin)
% PVL_EST_RS_PYSCH estimates Rs using multiple IV curves [1].
%
% Syntax
%   [Rs, errest] = pvl_est_Rs_Pysch(IVCurves, delI)
%   [Rs, errest] = pvl_est_Rs_Pysch(IVCurves, delI, Rsh, nNsVth, Io)
%
% Description
%   The method of Pysch [1] extends the Swanson method to use multiple 
%   IV curves. IV curves are assumed to be at different irradiance levels
%   but the same cell temperature. A point is selected on each IV curve where 
%   the current is less than Isc by the input delI. If optional arguments are
%   provided, the difference between the returned Rs value and the Rs
%   parameter for the single diode equation is estimated, see [2].
%
% Inputs:
%   IVCurves - structure array for IV curves including fields I and V
%   delI - offset from Isc to use
%   Rsh - (optional) a vector of Rsh in ohms for each IV curve
%   a - (optional) a vector of nNsVth for each IV curve
%   Io - (optional) a vector of dark current for each IV curve
%
% Output:
%   Rs - the series resistance value in ohms.
%   errest - the estimated difference between Rs and the series resistance
%     parameter for the single diode equation.
%
% References
%   [1] D. Pysch, A. Mette, S. W. Glunz, “A review and comparison of 
%   different methods to determine the series resistance of solar cells",
%   Solar. Energy Materials and Cells 91, pp. 1698-1706, 2007.
%
%   [2] C. Hansen and B. King, "Determining series resistance for
%   equivalent circuit models of a PV module", in 45th IEEE Photovoltaic
%   Specialist Conference, Waikoloa, HI, 2018.

p = inputParser;
addRequired(p, 'IVCurves', @isstruct);
addRequired(p, 'delI', @isnumeric);

addOptional(p, 'Rsh', NaN, @isnumeric);
addOptional(p, 'nNsVth', NaN, @isnumeric);
addOptional(p, 'Io', NaN, @isnumeric);

parse(p, IVCurves, delI, varargin{:});
IVCurves = p.Results.IVCurves;
delI = p.Results.delI;

Rsh = p.Results.Rsh;
nNsVth = p.Results.nNsVth;
Io = p.Results.Io;

% sort IV curves in order of increasing Isc

Isc = [IVCurves(:).Isc];
[~, Iscidx] = sort(Isc);

IVCurves = IVCurves(Iscidx);

pI = NaN(size(IVCurves)); pV = pI;
% estimate Rs
for k=1:length(IVCurves)
    % find current and voltage on each IV curve
    pI(k) = IVCurves(k).Isc - delI;
    if pI(k)>0
        pV(k) = interp1(IVCurves(k).I, IVCurves(k).V, pI(k));
    end
end
pI = pI(:); pV = pV(:);

X = [ones(size(pI(~isnan(pI)))) pI(~isnan(pI))];
Y = pV(~isnan(pI));
beta = X\Y;
Rs = -beta(2);

% % make Pysch Fig 4
% figure
% hold all
% for k=1:length(IVCurves)
%     plot(IVCurves(k).V, -IVCurves(k).I, '.')
%     plot(pV(k), -pI(k), 'bs')
% end
% plot(pV, -(1/beta(2)*pV - beta(1)/beta(2)), 'r-')

% for error estimate 
% create index for list of pairs of IV curves
[P, Q] = meshgrid(1:length(IVCurves), 1:length(IVCurves));
idx = [P(:), Q(:)];
% remove index pairs where i<=j
u = idx(:,1)<=idx(:,2);
idx = idx(~u, :);

% calculate weights
diffI = pI(idx(:,1)) - pI(idx(:,2));
w = diffI.^2 ./ nansum(nansum(diffI.^2));

if exist('a','var') && exist('Io','var') && exist('Rsh','var')
    
    v = Rsh(idx(:,1))>0 & Rsh(idx(:,2))>0 & Io>0 & nNsVth>0;

    % calculate estimated error
    delta = log(Rsh(idx(v,1))./Rsh(idx(v,2))) + ...
            log((delI.*Rsh(idx(v,2))./nNsVth + log(Io.*Rsh(idx(v,2))./nNsVth)) ./ ...
                (delI.*Rsh(idx(v,1))./nNsVth + log(Io.*Rsh(idx(v,1))./nNsVth)));
    delta2 = 1.04*( log(delI.*Rsh(idx(v,1)) + nNsVth.*log(Io.*Rsh(idx(v,1))./nNsVth)) ./ ...
                       (delI.*Rsh(idx(v,1)) + nNsVth.*log(Io.*Rsh(idx(v,1))./nNsVth))  - ...
                    log(delI.*Rsh(idx(v,2)) + nNsVth.*log(Io.*Rsh(idx(v,2))./nNsVth)) ./ ...
                       (delI.*Rsh(idx(v,2)) + nNsVth.*log(Io.*Rsh(idx(v,2))./nNsVth)) );
             
    errest = nansum(w(v).*nNsVth./diffI(v).*(delta + delta2));

else
    errest = NaN;
end


