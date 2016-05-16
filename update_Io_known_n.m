function [outIo] = update_Io_known_n(Rsh, Rs, nNsVth, Io, IL, Voc)

% UPDATE_IO_KNOWN_N adjusts Io to match Voc using other parameter values
%
% Syntax
%   [outIo] = update_Io_known_n(Rsh, Rs, nNsVth, Io, IL, Voc)
%
% Description
%   update_Io_known_n adjusts Io to match Voc using other parameter values,
%   i.e., Rsh (shunt resistance), Rs (series resistance), n (diode factor)
%   and IL (light current).  Io is updated iteratively 10 times or until
%   successive values are less than 0.000001% different.  The updating is
%   similar to Newton's method.
%
% Input:
%   Rsh - a vector of length N of values for shunt resistance (ohm)
%   Rs - a vector of length N of values for series resistance (ohm)
%   nNsVth - a vector of length N of values for the diode factor x thermal 
%     voltage for the module, equal to Ns (number of cells in series) 
%       x Vth (thermal voltage per cell).
%   Io - a vector of length N of initial values for Io (A)
%   IL - a vector of length N of values for light current IL (A)
%   Voc - a vector of length N of values for Voc (V)
%
% Output:
%   outIo - a vector of length N of updated values for Io
%
% Sources:
% [1] C. Hansen, Parameter Estimation for Single Diode Models of 
%     Photovoltaic Modules, Sandia National Laboratories Report SAND2015-XXXX
% [2] C. Hansen, Estimation of Parameters for Single Diode Models using
%     Measured IV Curves, Proc. of the 39th IEEE PVSC, June 2013.

eps = 1e-6;
niter = 10;
k = 1;
maxerr = 1;

tIo=Io(:);  % current estimate of Io

while maxerr>eps & k<niter
    % Predict Voc
    pVoc = pvl_V_from_I(Rsh,Rs,nNsVth,0,tIo,IL);
    % difference in Voc
    dVoc = pVoc - Voc;
    % Update Io
    next_Io = tIo.*(1+(2*dVoc)./(2*nNsVth - dVoc)); 
    % calculate maximum percent difference
    maxerr = max(abs(next_Io - tIo)./tIo)*100;
    tIo = next_Io;
    k = k+1;
end

outIo = tIo;

function [V] = pvl_V_from_I(Rsh, Rs, nNsVth, I, Io, Iphi)

% PVL_V_from_I Calculate voltage V from current I for a single diode equation 
%
% Syntax
%   V = pvl_V_from_I(Rsh, Rs, nNsVth, I, Io, Iphi)
% 
% Description
%   Calculate voltage V from current I using the single diode equation and 
%   parameter values Rsh, Rs, nNsVth, I, Io, Iphi.
%
% Inputs: 
%       Rsh - shunt resistance (ohm)
%       Rs - series resistance (ohm)
%       nNsVth - product of diode factor n, number of cells in series Ns
%           and cell thermal voltage Vth
%       I - current at which corresponding voltage will be computed (A)
%       Io - dark current (A)
%       Iphi - light current (A)
%
%   Inputs may be scalar or vectors, but all vectors must be of the same
%   length.
%
% Outputs:
%   V - voltage corresponding to I.
%
% Sources:
%
% [1] A. Jain, A. Kapoor, "Exact analytical solutions of the parameters of 
%     real solar cells using Lambert W-function", Solar Energy Materials 
%     and Solar Cells, 81 (2004) 269-277.
%
% See also PVL_LAMBERTW

% Generate the argument of the LambertW function
argW = (Io.*Rsh./nNsVth) .* ...
    exp(Rsh.*(-I+Iphi+Io)./(nNsVth));
inputterm = pvl_lambertw(argW); % Get the LambertW output
f = isnan(inputterm); % If argW is too big, the LambertW result will be NaN and we have to go to logspace

% Where f=NaN then the input argument (argW) is too large. It is necessary to go to logspace
if any(f)
    % Calculate the log(argW) if argW is really big
    logargW = log(Io) + log(Rsh) + Rsh.*(Iphi+Io-I)./(nNsVth)...
        - (log(nNsVth));
    
    % Three iterations of Newton-Raphson method to solve w+log(w)=logargW.
    % The initial guess is w=logargW. Where direct evaluation (above) results
    % in NaN from overflow, 3 iterations of Newton's method gives 
    % approximately 8 digits of precision.
    w = logargW;  
    for i=1:3
        w = w.*((1-log(w)+logargW)./(1+w));
    end;
    inputterm(f) = w(f);
end

% Eqn. 3 in Jain and Kapoor, 2004
V = -I.*(Rs + Rsh) + Iphi .* Rsh - nNsVth .* ...
    inputterm + Io .* Rsh;
