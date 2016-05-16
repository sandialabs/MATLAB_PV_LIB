function IAM = pvl_ashraeiam(b, theta)
% PVL_ASHRAEIAM Determine the incidence angle modifier using the ASHRAE
%   transmission model.
%
% Syntax
%   IAM = pvl_ashraeiam(b, theta)
%
% Description:
%   pvl_ashraeiam calculates the incidence angle modifier as developed in
%   [1], and adopted by ASHRAE (American Society of Heating, Refrigeration,
%   and Air Conditioning Engineers) [2]. The model has been used by model
%   programs such as PVSyst [3].
%
%   Note: For incident angles near 90 degrees, this model has a
%   discontinuity which has been addressed in this function.
%
% Input Parameters:
%   b - A parameter to adjust the modifier as a function of angle of
%     incidence. Typical values are on the order of 0.05 [3].
%   theta - The angle of incidence between the module normal vector and the
%     sun-beam vector in degrees. Theta must be a numeric scalar or vector.
%     For any values of theta where abs(theta)>90, IAM is set to 0. For any
%     values of theta where -90 < theta < 0, theta is set to abs(theta) and
%     evaluated. A warning will be generated if any(theta<0 or theta>90).
%     For values of theta near 90 degrees, the ASHRAE model may be above 1
%     or less than 0 due to the discontinuity of secant(theta). IAM values
%     outside of [0,1] are set to 0 and a warning is generated.
% 
% Output Parameters:
%   IAM - The incident angle modifier calculated as 1-b*(sec(theta)-1) as
%     described in [2,3]. IAM is a column vector with the same number of 
%     elements as the largest input vector.
%
% References:
%
% [1] Souka A.F., Safwat H.H., Determindation of the optimum orientations
%     for the double exposure flat-plate collector and its reflections,
%     Solar Energy vol .10, pp 170-174. 1966.
%
% [2] ASHRAE Standard 93-2010, Methods for Testing to Determine the
%     Thermal Performance of Solar Collectors
%
% [3] PVsyst Contextual Help. 
%     http://files.pvsyst.com/help/index.html?iam_loss.htm retrieved on
%     September 10, 2012
%
% See also 
%      PVL_GETAOI   PVL_EPHEMERIS   PVL_SPA   PVL_PHYSICALIAM


p = inputParser;
p.addRequired('bzero', @(x) (isnumeric(x) & all(x>=0) & isvector(x)));
p.addRequired('theta', @(x) (isnumeric(x) & isvector(x)));
p.parse(b, theta);

bzero = p.Results.bzero(:);
theta = p.Results.theta(:);

% theta(theta>-90 & theta<0) = abs(theta(theta>-90 & theta<0))
% First, any IAM(abs(theta)>=90) = 0
% Then, any IAM(IAM < 0 | IAM >1)=0 

if any(theta<0 | theta >= 90)
    % Alert the user of incorrect input angles and what we're going to do
    % about them
    warning(['Input incident angles <0 or >=90 detected in ',mfilename,'. ',...
        'For input angles with absolute value greater than 90, the ', ...
        'modifier is set to 0. For input angles between -90 and 0, the ',...
        'angle is changed to its absolute value and evaluated.']);
    % Set any negative input angles to be their absolute value
    theta(sign(theta)==-1)=abs(theta(sign(theta)==-1));
    

end

IAM = 1-bzero.*(secd(theta) - 1);
IAM(abs(theta)>90)=0;

if any((IAM > 1) | (IAM < 0))
    % Alert the user that we're truncating the output to remove the
    % discontinuity
    warning(['It seems that we have encountered a discontinuity in ',...
        mfilename,'. Any incident angle modifiers calculated to be less than 0 or '...
        'greather than 1 have been set to 0.']);
end
IAM((IAM > 1) | (IAM < 0))=0;
