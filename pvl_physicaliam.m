function IAM = pvl_physicaliam(K, L, n, theta)
% PVL_PHYSICALIAM Determine the incidence angle modifier using refractive 
%   index, glazing thickness, and extinction coefficient
%
% Syntax
%   IAM = pvl_physicaliam(K, L, n, theta)
%
% Description:
%   pvl_physicaliam calculates the incidence angle modifier as described in
%   De Soto et al. "Improvement and validation of a model for photovoltaic
%   array performance", section 3. The calculation is based upon a physical
%   model of absorbtion and transmission through a cover. Required
%   information includes, incident angle, cover extinction coefficient,
%   cover thickness
%
%   Note: The authors of this function believe that eqn. 14 in [1] is
%   incorrect. This function uses the following equation in its place:
%   theta_r = arcsin(1/n * sin(theta))
%
% Input Parameters:
%   K - The glazing extinction coefficient in units of 1/meters. Reference
%     [1] indicates that a value of  4 is reasonable for "water white"
%     glass. K must be a numeric scalar or vector with all values >=0. If K
%     is a vector, it must be the same size as all other input vectors.
%   L - The glazing thickness in units of meters. Reference [1] indicates
%     that 0.002 meters (2 mm) is reasonable for most glass-covered
%     PV panels. L must be a numeric scalar or vector with all values >=0. 
%     If L is a vector, it must be the same size as all other input vectors.
%   n - The effective index of refraction (unitless). Reference [1]
%     indicates that a value of 1.526 is acceptable for glass. n must be a 
%     numeric scalar or vector with all values >=0. If n is a vector, it 
%     must be the same size as all other input vectors.
%   theta - The angle of incidence between the module normal vector and the
%     sun-beam vector in degrees. Theta must be a numeric scalar or vector.
%     For any values of theta where abs(theta)>90, IAM is set to 0. For any
%     values of theta where -90 < theta < 0, theta is set to abs(theta) and
%     evaluated. A warning will be generated if any(theta<0 or theta>90).
% 
% Output Parameters:
%   IAM - The incident angle modifier as specified in eqns. 14-16 of [1].
%     IAM is a column vector with the same number of elements as the
%     largest input vector.
%
% References:
%
% [1] W. De Soto et al., "Improvement and validation of a model for
%     photovoltaic array performance", Solar Energy, vol 80, pp. 78-88,
%     2006.
%
% [2] Duffie, John A. & Beckman, William A.. ( © 2006). Solar Engineering 
%     of Thermal Processes, third edition. [Books24x7 version] Available 
%     from http://common.books24x7.com/toc.aspx?bookid=17160. 
%
% See also 
%      PVL_GETAOI   PVL_EPHEMERIS   PVL_SPA    PVL_ASHRAEIAM
%      PVL_MARTINRUIZIAM

% theta = incident angle in degrees
% n = effective refractive index
% K = glazing extinction coefficient (1/meters)
% L = glazing thickness (meters)
%
p=inputParser;
p.addRequired('K', @(x) (isnumeric(x) & all(x>=0) & isvector(x)));
p.addRequired('L', @(x) (isnumeric(x) & all(x>=0) & isvector(x)));
p.addRequired('n', @(x) (isnumeric(x) & all(x>=0) & isvector(x)));
p.addRequired('theta', @(x) (isnumeric(x) & isvector(x)));
p.parse(K, L, n, theta);

thetainput = p.Results.theta(:);
thetaused = thetainput;
K = p.Results.K(:);
n = p.Results.n(:);
L = p.Results.L(:);

% Check to see if any input angles are less than 0 or greater than 90.
if any(thetainput<0 | thetainput >90)
    % Alert the user of incorrect input angles and what we're going to do
    % about them
    warning(['Input incident angles <0 or >90 detected in ',mfilename,'. '...
        'For input angles with absolute value greater than 90, the ' ...
        'modifier is set to 0. For input angles between -90 and 0, the '...
        'angle is changed to its absolute value and evaluated.']);
    % Set any negative input angles to be their absolute value
    thetaused(sign(thetainput)==-1)=abs(thetainput(sign(thetainput)==-1));
    % Set any input angles with absolute value greater than 90 to 90. This
    % will be invalidated later, but it's a good step for right here.
    thetaused(abs(thetainput)>90)=90;
end

% Calculate theta_r using what we believe to be a corrected version of
% equation 14 in [1]
thetar_deg = asind(1./n.*sind(thetaused));

% Calculate equation 15 to get Tau
tau = exp(-1.*(K.*L./cosd(thetar_deg))).*...
    (1 - 0.5 .* ( ((sind(thetar_deg-thetaused)).^2)./((sind(thetar_deg+thetaused)).^2) + ...
    ((tand(thetar_deg-thetaused)).^2)./((tand(thetar_deg+thetaused)).^2)));

zeroang = 0.000001; % Pick an angle very close to 0, but not 0.
thetar_deg0 = asind(1./n.*sind(zeroang)); % Find theta_r(0) using the small angle

% Find Tau(0)
tau0 = exp(-1.*(K.*L./cosd(thetar_deg0))).*...
    (1 - 0.5 .* ( ((sind(thetar_deg0-zeroang)).^2)./((sind(thetar_deg0+zeroang)).^2) + ...
    ((tand(thetar_deg0-zeroang)).^2)./((tand(thetar_deg0+zeroang)).^2)));

IAM = tau./tau0; % Equation 16
IAM(thetainput==0)=1; % Normal incidence should have an incident angle modifier of 1
IAM((abs(thetainput)>90) | (IAM < 0))=0; % Set the modifier to 0 for any input angles with absolute value > 90
end