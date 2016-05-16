function [outRsh] = update_Rsh_fixed_pt(Rsh, Rs, Io, IL, nNsVth, Imp, Vmp)

% UPDATE_RSH_FIXED_PT adjusts Rsh to match Vmp using other parameter values
%
% Syntax
%   [outRsh] = update_Rsh_fixed_pt(Rsh, Rs, Io, IL, nNsVth, Imp, Vmp)
%
% Description
%   update_Rsh_fixed_pt adjusts Rsh to match Vmp using other parameter values,
%   i.e., Rs (series resistance), n (diode factor), Io (dark current), 
%   and IL (light current).  Rsh is updated iteratively using a fixed point
%   expression obtained from combining Vmp=Vmp(Imp) (using the analytic
%   solution to the single diode equation) and dP/dI = 0 at Imp.  500
%   iterations are performed because convergence can be very slow.
%
% Input:
%   Rsh - a vector of length N of initial values for shunt resistance (ohm)
%   Rs - a vector of length N of values for series resistance (ohm)
%   Io - a vector of length N of values for Io (A)
%   IL - a vector of length N of values for light current IL (A)
%   nNsVth - a vector of length N of values for the diode factor x thermal 
%     voltage for the module, equal to Ns (number of cells in series) 
%       x Vth (thermal voltage per cell).
%   Imp - a vector of length N of values for Imp (V)
%   Vmp - a vector of length N of values for Vmp (V)
%
% Output:
%   outRsh - a vector of length N of updated values for Rsh
%
% Sources:
% [1] C. Hansen, Parameter Estimation for Single Diode Models of 
%     Photovoltaic Modules, Sandia National Laboratories Report SAND2015-XXXX

niter = 500;
x1=Rsh(:);

for k=1:niter
    [~, z]=calc_theta_phi_exact(Imp,IL,Vmp,Io,nNsVth,Rs,x1);
    next_x1 = (1+z)./z.*((IL+Io).*x1./Imp-nNsVth.*z./Imp-2*Vmp./Imp);
    x1 = next_x1;
end

outRsh = x1;
