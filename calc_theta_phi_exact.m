function [theta phi] = calc_theta_phi_exact(Imp, IL, Vmp, Io, nNsVth, Rs, Rsh)

% CALC_THETA_PHI_EXACT computes Lambert W values appearing in the analytic
%   solutions to the single diode equation for the max power point.
%
% Syntax
%   [theta phi] = calc_theta_phi_exact(Imp, IL, Vmp, Io, nNsVth, Rs, Rsh)
%
% Description
%   calc_theta_phi_exact calculates values for the Lambert W function 
%   which are used in the analytic solutions for the single diode equation
%   at the maximum power point.
%   For V=V(I), phi = W(Io*Rsh/n*Vth * exp((IL + Io - Imp)*Rsh/n*Vth))
%   For I=I(V), theta = W(Rs*Io/n*Vth * Rsh/(Rsh+Rs) 
%                             * exp(Rsh/(Rsh+Rs)*((Rs(IL+Io) + V)/n*Vth))
%
% Input:
%   Imp - a vector of length N of values for Imp (A)
%   IL - a vector of length N of values for light current IL (A)
%   Vmp - a vector of length N of values for Vmp (V)
%   Io - a vector of length N of values for Io (A)
%   nNsVth - a vector of length N of values for the diode factor x thermal 
%     voltage for the module, equal to Ns (number of cells in series) 
%       x Vth (thermal voltage per cell).
%   Rs - a vector of length N of values for series resistance (ohm)
%   Rsh - a vector of length N of values for shunt resistance (ohm)
%
% Output:
%   theta - a vector of values for the Lambert W function for solving
%      I=I(V)
%   phi - a vector of values for the Lambert W function for solving V=V(I)
%
% Sources:
%   [1] C. Hansen, Parameter Estimation for Single Diode Models of 
%   Photovoltaic Modules, Sandia National Laboratories Report SAND2015-XXXX
%   [2] A. Jain, A. Kapoor, "Exact analytical solutions of the parameters of 
%     real solar cells using Lambert W-function", Solar Energy Materials 
%     and Solar Cells, 81 (2004) 269-277.

% argument for Lambert W function involved in V = V(I)
% [1] Eq. 12; [2] Eq. 3
argw = Rsh.*Io./nNsVth .*exp(Rsh.*(IL + Io - Imp)./nNsVth);
u = argw>0;
w(~u)=NaN;
tmp = pvl_lambertw(argw(u));
ff = isnan(tmp);

% NaN where argw overflows. Switch to log space to evaluate
if any(ff)
    logargW = log(Rsh(u)) + log(Io(u)) - log(nNsVth(u)) + Rsh(u).*(IL(u) + Io(u) - Imp(u))./nNsVth(u);
    % Three iterations of Newton-Raphson method to solve w+log(w)=logargW.
    % The initial guess is w=logargW. Where direct evaluation (above) results
    % in NaN from overflow, 3 iterations of Newton's method gives 
    % approximately 8 digits of precision.
    x = logargW;  
    for i=1:3
        x = x.*((1-log(x)+logargW)./(1+x));
    end;
    tmp(ff) = x(ff);
end
w(u) = tmp;
phi = w(:);

% argument for Lambert W function involved in I=I(V)
% [1] Eq. 11; [2] Eq. 2
argw = Rsh./(Rsh+Rs).*Rs.*Io./nNsVth .*exp(Rsh./(Rsh+Rs).*(Rs.*(IL+Io)+Vmp)./nNsVth);
u = argw>0;
w(~u)=NaN;
tmp = pvl_lambertw(argw(u));
ff = isnan(tmp);

% NaN where argw overflows. Switch to log space to evaluate
if any(ff)
    logargW = log(Rsh(u)./(Rsh(u)+Rs(u))) + log(Rs(u)) + log(Io(u)) - log(nNsVth(u)) + (Rsh(u)./(Rsh(u)+Rs(u))).*(Rs(u).*(IL(u)+Io(u)) + Vmp(u))./nNsVth(u);
    % Three iterations of Newton-Raphson method to solve w+log(w)=logargW.
    % The initial guess is w=logargW. Where direct evaluation (above) results
    % in NaN from overflow, 3 iterations of Newton's method gives 
    % approximately 8 digits of precision.
    x = logargW;  
    for i=1:3
        x = x.*((1-log(x)+logargW)./(1+x));
    end;
    tmp(ff) = x(ff);
end
w(u) = tmp;
theta = w(:);
