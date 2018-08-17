function [IL, I0, Rsh, Rs] = pvl_est_diode_params_simple(IVCurve, nNsVth)
% PVL_EST_SINGLE_DIODE_PARAM_SIMPLE fits the single diode equation to data
% for a single IV curve.
%
% Syntax
%   [IL, I0, Rsh, Rs] = pvl_est_single_diode_param_simple(IVCurve, nNsVth)
%
% Description
%   pvl_est_single_diode_param_simple uses a sequential technique described
%   in [1] to fit the single diode equation to data for a single IV curve.
%   The method here is a simplification of that found in [2] and coded in
%   est_single_diode_param. The diode factor n is not determined in this
%   function. A value for n can be found by regression between Voc and 
%   log(Ee) for a range of effective irradiance Ee, as coded in e.g. 
%   pvl_desoto_parameter_estimation.
%
% Input:
%   IVCurve - A structure with the following fields:
%      IVCurve.Isc - short circuit current in amperes.
%      IVCurve.Voc - open circuit voltage in volts.
%      IVCurve.Imp - current at maximum power point in amperes. 
%      IVCurve.Vmp - voltage at maximum power point in volts.
%      IVCurve.Pmp - power at maximum power point in watts.
%      IVCurve.V - vector of voltage in volts. 
%      IVCurve.I - vector of current in amperes.
%   nNsVth - the product n (diode factor) x Ns (cells in series)
%     x Vth (thermal voltage per cell) for both IV curves.
%
% Output:
%   IL - the light current (A) for the IV curve
%   I0 - the dark current (A) for the IV curve
%   Rsh - shunt resistance (ohm) for the IV curve
%   Rs - series resistance (ohm) for the IV curve
%
% References
%   [1] C. Hansen and B. King, "Determining series resistance for
%   equivalent circuit models of a PV module", in 45th IEEE Photovoltaic
%   Specialist Conference, Waikoloa, HI, 2018.
%
%   [2] C. Hansen, Parameter Estimation for Single Diode Models of 
%   Photovoltaic Modules, Sandia National Laboratories Report SAND2015-2065

I = IVCurve.I(:);
V = IVCurve.V(:);

% rule of thumb for portion of IV curve I ~ IL + I0 - (V + I * Rs)/ Rsh
% where exponential term has small influence. Isc stands in here for IL.

Vlim = nNsVth*log(5*10^2) - IVCurve.Isc;

idx = find(V - Vlim > 0, 1);
tRsh = -1;
% a bit of protection against IV curves where measured current is slightly
% greater than Isc at some voltages
while tRsh<0 && idx<length(V)
    X = V(1:idx);
    Y = IVCurve.Isc - I(1:idx);
    beta = [ones(size(X)) X]\Y;
    tRsh = 1/beta(2);
    idx = idx + 5;
end

% initial guess at IL
IL = IVCurve.Isc;

for k=1:5
    % calculate I0
    if IVCurve.Voc/nNsVth < (log(realmax)-3)
        I0 = (IL - IVCurve.Voc/tRsh)/(exp(IVCurve.Voc/nNsVth) - 1);
    else
        logI0 = log(IL - Voc/tRsh) - IVCurve.Voc/nNsVth;
        I0 = exp(logI0);
    end

    % using values for nNsVth, Rsh, IL and I0, calculate Rs at a point
    % midway between Vmp and Voc
    idx = find(V - (IVCurve.Voc + IVCurve.Vmp)/2 > 0, 1);
    tV = V(idx);
    tI = I(idx);
    W = calc_phi_exact(tI, IL, I0, nNsVth, tRsh);
    Rs = ((IL + I0 - tI)*tRsh - tV - nNsVth*W)/tI;

    % update IL
    IL = IVCurve.Isc*(1 + Rs/tRsh);
end

Rsh = tRsh;

end


function [W] = calc_phi_exact(I, IL, Io, a, Rsh)

% calculates W(phi) where phi is the argument of the
% Lambert W function in V = V(I) at I=Imp ([2], Eq. 3).  Formula for
% phi is given in code below as argw.

% phi
argw = Rsh.*Io./a .*exp(Rsh.*(IL + Io - I)./a);

% Screen out any negative values for argw
u = argw>0;
W(~u)=NaN;

tmp = pvl_lambertw(argw(u));

ff = isnan(tmp);

% take care of any numerical overflow by evaluating log(W(phi))
if any(ff)
    logargW = log(Rsh(u)) + log(Io(u)) - log(a(u)) + Rsh(u).*(IL(u) + Io(u) - I(u))./a(u);
    % Three iterations of Newton-Raphson method to solve w+log(w)=logargW.
    % The initial guess is w=logargW. Where direct evaluation (above) results
    % in NaN from overflow, 3 iterations of Newton's method gives 
    % approximately 8 digits of precision.
    x = logargW;  
    for i=1:5
        x = x.*((1-log(x)+logargW)./(1+x));
    end
    tmp(ff) = x(ff);
end

W(u) = tmp;

end


    
    
