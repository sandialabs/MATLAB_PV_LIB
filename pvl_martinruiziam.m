function IAM = pvl_martinruiziam(ar, theta)
% PVL_PHYSICALIAM Determine the incidence angle modifier using the Martin 
%   and Ruiz incident angle model
%
% Syntax
%   IAM = pvl_martinruiziam(ar, theta)
%
% Description:
%   pvl_martinruiziam calculates the incidence angle modifier (angular
%   factor) as described by Martin and Ruiz in [1]. The information
%   required is the incident angle (theta) and the angular losses
%   coefficient (ar). Please note that [1] has a corrigendum which makes the
%   document much simpler to understand.
%
%
% Input Parameters:
%   ar - The angular losses coefficient described in equation 3 of [1].
%     This is an empirical dimensionless parameter. Values of ar are
%     generally on the order of 0.08 to 0.25 for flat-plate PV modules. ar
%     must be a numeric scalar or vector with all values > 0. If ar
%     is a vector, it must be the same size as all other input vectors.
%   theta - The angle of incidence between the module normal vector and the
%     sun-beam vector in degrees. Theta must be a numeric scalar or vector.
%     For any values of theta where abs(theta)>90, IAM is set to 0. For any
%     values of theta where -90 < theta < 0, theta is set to abs(theta) and
%     evaluated. A warning will be generated if any(theta<0 or theta>90).
% 
% Output Parameters:
%   IAM - The incident angle modifier from [1]. The incident angle modifier
%     is defined as [1-exp(-cos(theta/ar))] / [1-exp(-1/ar)], which is
%     presented as AL(alpha) = 1 - IAM in equation 4 of [1]. Thus IAM is
%     equal to 1 at theta = 0, and equal to 0 at theta = 90. IAM is a
%     column vector with the same number of elements as the largest input
%     vector.
%
% References:
%
% [1] N. Martin and J. M. Ruiz, "Calculation of the PV modules angular
%     losses under field conditions by means of an analytical model", Solar
%     Energy Materials & Solar Cells, vol. 70, pp. 25-38, 2001.
%
% [2] N. Martin and J. M. Ruiz, "Corrigendum to 'Calculation of the PV
%     modules angular losses under field conditions by means of an
%     analytical model'", Solar Energy Materials & Solar Cells, vol. 110,
%     pp. 154, 2013.
%
% See also 
%      PVL_GETAOI   PVL_EPHEMERIS   PVL_SPA    PVL_ASHRAEIAM 
%      PVL_PHYSICALIAM

% theta = incident angle in degrees
% ar = angular losses factor
%
p=inputParser;
p.addRequired('ar', @(x) (isnumeric(x) & all(x>0) & isvector(x)));
p.addRequired('theta', @(x) (isnumeric(x) & isvector(x)));
p.parse(ar, theta);

thetainput = p.Results.theta(:);
thetaused = thetainput;
ar = p.Results.ar(:);


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

% Calculate the incident angle modifier using equation 4 in [1]. Note that
% our definition of incident angle modifier (IAM) is:
% IAM = -1 * AL(alpha) + 1

IAM = (1-exp(-cosd(thetaused)./ar)) ./ (1-exp(-1/ar));

IAM((abs(thetainput)>90) | (IAM < 0))=0; % Set the modifier to 0 for any input angles with absolute value > 90
end