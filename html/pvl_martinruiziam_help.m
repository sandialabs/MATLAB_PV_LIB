%% pvl_martinruiziam 
% Determine the incidence angle modifier using the Martin and Ruiz model.
%
%% Syntax
% * |IAM = pvl_martinruiziam(ar, theta)|
%
%% Description:
% |pvl_martinruiziam| calculates the incidence angle modifier (IAM) 
% as described in [1] and clarified in [2].  The incidence angle
% modifier quantifies the fraction of direct irradiance incidence on a
% module's surface that is not reflected away.
% The incidence angle modifier is computed as 
%
% $$IAM = {\frac{1 - \exp(-\cos(\theta /a_r))} {1 - \exp(-1/a_r)}} $$,
%
% Thus $IAM = 1$ at $\theta = 0$, and $IAM = 0$ at $\theta = 90$.
%
%
%% Inputs
%%
% * *|ar|* - The angular losses coefficient described in equation 3 of [1].
%     |ar| is an empirical dimensionless parameter with values 
%     generally on the order of 0.08 to 0.25 for flat-plate PV modules.
%     |ar| must be a numeric scalar or vector with all values > 0. If |ar|
%     is a vector, it must be the same size as all other input vectors.
% * *|theta|* - The angle of incidence between the module normal vector and the
%     sun vector in degrees. |theta| must be a numeric scalar or vector.
%     For any values of |theta| where |abs(theta)>90|, |IAM| is set to 0. For any
%     values of |theta| where -90 < |theta| < 0, |theta| is set to |abs(theta)| and
%     evaluated. A warning will be generated if |any(theta<0 or theta>90)|.
% 
%% Output
% * *|IAM|* - The incident angle modifier.  |IAM| is a
%     column vector with the same number of elements as the largest input
%     vector.
%
%% Example
% This example plots the IAM for glass over a range of incident angles.
ar = 0.15;       % empirical losses coefficient of 0.15
theta = 0:90;    %incident angle in degrees
IAM = pvl_martinruiziam(ar, theta);
figure
plot(theta,IAM)
xlabel('Incident Angle (deg)')
ylabel('IAM')
title('Martin and Ruiz IAM Model Example')
%% References
%
% [1] N. Martin and J. M. Ruiz, 2001. Calculation of the PV modules angular
%     losses under field conditions by means of an analytical model, Solar
%     Energy Materials & Solar Cells, vol. 70, pp. 25-38.
%
% [2] N. Martin and J. M. Ruiz, 2013. Corrigendum to 'Calculation of the PV
%     modules angular losses under field conditions by means of an
%     analytical model', Solar Energy Materials & Solar Cells, vol. 110,
%     pp. 154.
%
%% See also 
% <pvl_getaoi_help.html |pvl_getaoi|> ,           
% <pvl_ephemeris_help.html |pvl_ephemeris|> ,
% <pvl_spa_help.html |pvl_spa|> ,
% <pvl_ashraeiam_help.html |pvl_ashraeiam|> ,
% <pvl_physicaliam_help.html |pvl_physicaliam|> 

%%
% Copyright 2014 Sandia National Laboratories
