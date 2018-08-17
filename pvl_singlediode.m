function [Result] = pvl_singlediode(IL, I0, Rs, Rsh, nNsVth, varargin)
% PVL_SINGLEDIODE Solve the single-diode model to obtain a photovoltaic IV curve
%
% Syntax
%   [Result] = pvl_singlediode(IL, I0, Rs, Rsh, nNsVth)
%   [Result] = pvl_singlediode(IL, I0, Rs, Rsh, nNsVth, NumPoints)
%   [Result] = pvl_singlediode(IL, I0, Rs, Rsh, nNsVth, Vpoints)
%
% Description
%   pvl_singlediode solves the single diode equation [1]:
%   I = IL - I0*[exp((V+I*Rs)/(nNsVth))-1] - (V + I*Rs)/Rsh
%   for I and V when given IL, I0, Rs, Rsh, and nNsVth (nNsVth = n*Ns*Vth) which
%   are described later. pvl_singlediode returns a struct which contains
%   the 5 points on the I-V curve specified in SAND2004-3535 [3], and can
%   optionally provide a full IV curve with a user-defined number of
%   points. If all IL, I0, Rs, Rsh, and nNsVth are scalar, a single curve
%   will be returned, if any are vectors (of the same length), multiple IV
%   curves will be calculated.
%
% Input Parameters:
%   IL - Light-generated current (photocurrent) in amperes under desired IV
%      curve conditions. May be a scalar or vector, but vectors must be of
%      same length as all other input vectors.
%   I0 - Diode saturation current in amperes under desired IV curve
%      conditions. May be a scalar or vector, but vectors must be of same
%      length as all other input vectors.
%   Rs - Series resistance in ohms under desired IV curve conditions. May
%      be a scalar or vector, but vectors must be of same length as all
%      other input vectors.
%   Rsh - Shunt resistance in ohms under desired IV curve conditions. May
%      be a scalar or vector, but vectors must be of same length as all
%      other input vectors.
%   nNsVth - the product of three components. 1) The usual diode ideal 
%      factor (n), 2) the number of cells in series (Ns), and 3) the cell 
%      thermal voltage under the desired IV curve conditions (Vth).
%      The thermal voltage of the cell (in volts) may be calculated as 
%      k*Tcell/q, where k is Boltzmann's constant (J/K), Tcell is the
%      temperature of the p-n junction in Kelvin, and q is the elementary 
%      charge of an electron (coulombs). nNsVth may be a scalar or vector, 
%      but vectors must be of same length as all other input vectors.
%   NumPoints - Number of points to compute in each IV curve (optional).
%      Must be a finite scalar value. Non-integer values will be rounded 
%      to the next highest integer (ceil). The default
%      value is 0, resulting in no calculation of IV points other than
%      those specified in [3].
%   Vpoints - number of points, or values of voltage at which to return IV curve points
%      (optional). Must be a scalar, 1d array or a cell array of the same length 
%       as all other input vectors. If a scalar, then the specified number of points
%       is computed for each IV curve.  If a 1d array then points in that array
%       between 0 and Voc are returned for each IV curve.  If a cell array then
%       each element of the cell array contains a 1d array of voltage values
%       at which points are returned.
%
% Output:
%   Result - A structure array with the following fields. All fields have the
%   same number of rows as the largest input vector.
%      Result.Isc - Column vector of short circuit current in amperes.
%      Result.Voc - Column vector of open circuit voltage in volts.
%      Result.Imp - Column vector of current at maximum power point in amperes. 
%      Result.Vmp - Column vector of voltage at maximum power point in volts.
%      Result.Pmp - Column vector of power at maximum power point in watts.
%      Result.Ix - Column vector of current, in amperes, at V = 0.5*Voc.
%      Result.Ixx - Column vector of current, in amperes, at V = 0.5*(Voc+Vmp).
%      Result.V - Cell array of voltages in volts. Element n is an array of 
%        voltages running from V=0 to V=Voc. Voltages are either equally
%        spaced if NumPoints is input, or at values specified in input
%        Vpoints.
%      Result.I - Cell array of currents in amperes. Element n is an array
%        of current values corresponding to the voltage values in element n
%        of Result.V.
%
% Notes:
%    1) The solution employed to solve the implicit diode equation utilizes
%       the Lambert W function to obtain an explicit function of V=f(i) and
%       I=f(V) as shown in [2].
%
% Sources:
%
% [1] S.R. Wenham, M.A. Green, M.E. Watt, "Applied Photovoltaics" 
%     ISBN 0 86758 909 4
%
% [2] A. Jain, A. Kapoor, "Exact analytical solutions of the parameters of 
%     real solar cells using Lambert W-function", Solar Energy Materials 
%     and Solar Cells, 81 (2004) 269-277.
%
% [3] D. King et al, "Sandia Photovoltaic Array Performance Model",
%     SAND2004-3535, Sandia National Laboratories, Albuquerque, NM
%
% See also
%   PVL_CALCPARAMS_DESOTO  PVL_CALCPARAMS_PVSYST
%

eps = 1e-6; % tolerance (fraction) for deciding two voltage values are the same

% Vth is the thermal voltage of the cell in VOLTS
% n is the usual diode ideality factor
% nNsVth is n*Ns*Vth

p = inputParser;
p.addRequired('IL',@(x) all(x>=0) & isnumeric(x) & isvector(x) );
p.addRequired('I0', @(x) all(x>=0) & isnumeric(x) & isvector(x));
p.addRequired('Rs', @(x) all(x>=0) & isnumeric(x) & isvector(x));
p.addRequired('Rsh', @(x) all(x>=0) & isnumeric(x) & isvector(x));
p.addRequired('nNsVth',@(x) all(x>=0) & isnumeric(x) & isvector(x) );
p.addOptional('Points', NaN, @(x) (isscalar(x) && isfinite(x)) || ...
                (isvector(x) && all(isfinite(x))) || ...
                (iscell(x)));
p.parse(IL, I0, Rs, Rsh, nNsVth, varargin{:});

% Make all inputs into column vectors
IL=p.Results.IL(:);
I0=p.Results.I0(:);
Rs = p.Results.Rs(:);
Rsh = p.Results.Rsh(:);
nNsVth = p.Results.nNsVth(:);
Points = p.Results.Points;

% Ensure that all input values are either 1) scalars or 2) vectors of equal
% length. The "vector" part of the check is performed by the input parser,
% this merely checks that they are of equal length.
VectorSizes = [numel(IL), numel(I0), numel(Rs), numel(Rsh), numel(nNsVth)];
MaxVectorSize = max(VectorSizes);
if not(all((VectorSizes==MaxVectorSize) | (VectorSizes==1)))
    error('Input vectors IL, I0, Rs, Rsh, and nNsVth must either be scalars or vectors of the same length.');
end

% Ensure that if a cell array of voltage values is supplied, that one
% element is provided for each IV curve

if iscell(Points) && not(numel(Points)==MaxVectorSize)
    error('A cell array of voltage was input.  The length of the cell array must equal the length of all other vectors')
end

% If any input variable is not a scalar, then make any scalar input values
% into a column vector of the correct size.
if (MaxVectorSize > 1 && any(VectorSizes == 1))
    IL = IL.*ones(MaxVectorSize , 1);
    I0 = I0.*ones(MaxVectorSize , 1);
    Rs = Rs.*ones(MaxVectorSize , 1);
    Rsh = Rsh.*ones(MaxVectorSize , 1);
    nNsVth = nNsVth.*ones(MaxVectorSize , 1);
end

Imax = zeros(MaxVectorSize, 1);
Pmp = zeros(MaxVectorSize, 1);
Vmax = zeros(MaxVectorSize, 1);
Ix = zeros(MaxVectorSize, 1);
Ixx = zeros(MaxVectorSize, 1);
Voc = zeros(MaxVectorSize, 1);
Isc = zeros(MaxVectorSize, 1);

u = IL>0;

% Find Isc using Lambert W
Isc(u) = I_from_V(Rsh(u), Rs(u), nNsVth(u), 0, I0(u), IL(u));

% Find Voc using Lambert W
Voc(u) = V_from_I(Rsh(u), Rs(u), nNsVth(u), 0, I0(u), IL(u));


% Calculate I, V and P at the maximum power point

[Imax(u), Vmax(u), Pmp(u)]=calc_Pmp_bisect(IL(u),I0(u),nNsVth(u),Rs(u),Rsh(u));


% Find Ix and Ixx using Lambert W
Ix(u) = I_from_V(Rsh(u), Rs(u), nNsVth(u), 0.5*Voc(u), I0(u), IL(u));
Ixx(u) = I_from_V(Rsh(u), Rs(u), nNsVth(u), 0.5*(Voc(u)+Vmax(u)), I0(u), IL(u));

% Wrap up the results into the Result struct
Result.Voc = Voc;
Result.Vmp = Vmax;
Result.Imp = Imax;
Result.Ix = Ix;
Result.Ixx = Ixx;
Result.Pmp = Pmp;
Result.Isc = Isc;


if ~isnan(Points)
    % IV curve points are requested
    Result.I = cell(MaxVectorSize,1);
    Vpoints = cell(MaxVectorSize,1);

    for i=1:MaxVectorSize
        if u(i)
            if isscalar(Points)
                % Take care of any pesky non-integers by rounding up
                % User specified NumPoints rather than Vpoints. Create
                % Vpoints, a cell array with an array of voltage values for each IV curve.
                NumPoints = ceil(Points);
                Vpoints(i) = { (Voc(i))*(0:1/(NumPoints-1):1) }; 
            elseif isvector(Points)
                % Same points for each curve
                epsV = eps*Voc(i);
                Vpoints(i) = { Points(Points>=0 & Points<=Voc(i)) };
            elseif iscell(Points)
                Vpoints(i) = Points(i);
            end
            
            tmpI = I_from_V(Rsh(i), Rs(i), nNsVth(i), Vpoints{i}, I0(i), IL(i));
            
        else
            Vpoints(i) = { zeros(NumPoints) };
            tmpI = [];
        end
        Result.I(i) = { tmpI };

    end

    Result.V = cell2mat(Vpoints);
    Result.I = cell2mat(Result.I);
end


end


function [V]=V_from_I(Rsh, Rs, nNsVth, I, Io, Iphi)

% calculates V from I per Eq 3 Jain and Kapoor 2004
% uses Lambert W implemented in pvl_lambertw.m
% Rsh, nVth, I, Io, Iphi can all be vectors
% Rs can be a vector, but should be a scalar

% Generate the argument of the LambertW function
argW = (Io.*Rsh./nNsVth) .* ...
    exp(Rsh.*(-I+Iphi+Io)./(nNsVth));
inputterm = pvl_lambertw(argW); % Get the LambertW output
f = isnan(inputterm); % If argW is too big, the LambertW result will be NaN and we have to go to logspace

% If it is necessary to go to logspace
if any(f)
    % Calculate the log(argW) if argW is really big
    logargW = log(Io) + log(Rsh) + Rsh.*(Iphi+Io-I)./(nNsVth)...
        - (log(nNsVth));
    
    % Three iterations of Newton-Raphson method to solve w+log(w)=logargW.
    % The initial guess is w=logargW. Where direct evaluation (above) results
    % in NaN from overflow, 3 iterations of Newton's method gives 
    % approximately 8 digits of precision.
    w = logargW;
    K = round(log10(w));
    for i=1:(3*K)
        w = w.*((1-log(w)+logargW)./(1+w));
    end;
    inputterm(f) = w(f);
end

% Eqn. 3 in Jain and Kapoor, 2004. modified to better handle very large Rsh values

V = Rsh.*(-I.*(Rs./Rsh + 1) + Iphi - nNsVth./Rsh .* ...
    inputterm + Io);

end

function [I]=I_from_V(Rsh, Rs, nNsVth, V, Io, Iphi)

% calculates I from V per Eq 2 Jain and Kapoor 2004
% uses Lambert W implemented in pvl_lambertw.m
% Rsh, nVth, V, Io, Iphi can all be vectors
% Rs can be a vector, but should be a scalar

argW = Rs.*Io.*Rsh.*...
    exp(Rsh.*(Rs.*(Iphi+Io)+V)./(nNsVth.*(Rs+Rsh)))./...
    (nNsVth.*(Rs + Rsh));
inputterm = pvl_lambertw(argW);

% Eqn. 4 in Jain and Kapoor, 2004
I = -V./(Rs + Rsh) - (nNsVth./Rs) .* inputterm + ...
    Rsh.*(Iphi + Io)./(Rs + Rsh);
end

function [W] = calc_phi_exact(Imp, IL, Io, a, Rsh)

% calculates W(phi) where phi is the argument of the
% Lambert W function in V = V(I) at I=Imp ([2], Eq. 3).  Formula for
% phi is given in code below as argw.

% phi
argw = Rsh.*Io./a .*exp(Rsh.*(IL + Io - Imp)./a);

% Screen out any negative values for argw
u = argw>0;
W(~u)=NaN;

tmp = pvl_lambertw(argw(u));

ff = isnan(tmp);

% take care of any numerical overflow by evaluating log(W(phi))
if any(ff)
    logargW = log(Rsh(u)) + log(Io(u)) - log(a(u)) + Rsh(u).*(IL(u) + Io(u) - Imp(u))./a(u);
    % Three iterations of Newton-Raphson method to solve w+log(w)=logargW.
    % The initial guess is w=logargW. Where direct evaluation (above) results
    % in NaN from overflow, 3 iterations of Newton's method gives 
    % approximately 8 digits of precision.
    x = logargW; 
    K = round(log10(x));
    for i=1:(3*K)
        x = x.*((1-log(x)+logargW)./(1+x));
    end;
    tmp(ff) = x(ff);
end;

W(u) = tmp;

end

function Imp=calc_Imp_bisect(Iph,Io,a,Rs,Rsh)

% calculates the value of Imp (current at maximum power point) for an IV
% curve with parameters Iph, Io, a, Rs, Rsh.  Imp is found as the value of
% I for which g(I)=dP/dV (I) = 0.

% Set up lower and upper bounds on Imp
A = 0*Iph;
B = Iph+Io;

% Detect when lower and upper bounds are not consistent with finding the
% zero of dP/dV

gA=g(A,Iph,Io,a,Rs,Rsh);
gB=g(B,Iph,Io,a,Rs,Rsh);

if any(gA.*gB>0)
    % where gA*gB>0, then there is a problem with the IV curve parameters.
    % In the event where gA and gB have the same sign, alert the user with
    % a warning and replace erroneous cases with NaN
    errorvalues = gA .* gB > 0;
    warning(['Warning: pvl_singlediode has found at least one case where' ...
        ' the single diode parameters are such that dP/dV may not have' ...
        ' a zero. A NaN value has been reported for all such cases.'])
    A(errorvalues) = NaN; % This will set Imp values where gA*gB>0 to NaN
end

% midpoint is initial guess for Imp
p = (A+B)./2;
err = g(p,Iph,Io,a,Rs,Rsh);  % value of dP/dV at initial guess p

while max(abs(B-A))>1e-6   % set precision of estimate of Imp to 1e-6 (A)
    gA=g(A,Iph,Io,a,Rs,Rsh); % value of dP/dV at left endpoint
    u=(gA.*err)<0;
    B(u)=p(u);
    A(~u)=p(~u);
    p=(A+B)/2;
    err = g(p,Iph,Io,a,Rs,Rsh);
end;
Imp = p;

end

function y=g(I,Iph,Io,a,Rs,Rsh)

% calculates dP/dV exactly, using p=I*V=I*V(I), where V=V(I) uses the
% Lambert's W function W(phi) ([2], Eq. 3).

[z]=calc_phi_exact(I,Iph,Io,a,Rsh);  % calculate W(phi)
z = z(:);

% calculate dP/dV
y = (Iph+Io-2*I).*Rsh - 2*I.*Rs - a.*z +I.*Rsh.*z./(1+z);

end

function [Imp, Vmp, Pmp]=calc_Pmp_bisect(Iph,Io,a,Rs,Rsh)

% Returns Imp, Vmp, Pmp for the IV curve described by input parameters.
% Vectorized.

Imp=calc_Imp_bisect(Iph,Io,a,Rs,Rsh);  % find Imp
[z]=calc_phi_exact(Imp,Iph,Io,a,Rsh);  % Calculate W(phi) at Imp, where W is Lambert's W function and phi is its argument ([2], Eq. 3)
z = z(:);
Vmp=(Iph+Io-Imp).*Rsh - Imp.*Rs - a.*z;  % Compute V from Imp and W(phi)
Pmp = Vmp.*Imp;
end