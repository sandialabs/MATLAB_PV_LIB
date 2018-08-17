function [Rloss_Beam,Rloss_Iso,Rloss_Albedo] = pvl_iam_martinruiz_components(SurfTilt, AOI, varargin)

% PVL_IAM_MARTINRUIZ_COMPONENTS calculates the reflection loss for direct beam,
% isotropic diffuse, and ground-reflected albedo light.
%
% Syntax
%   |[Rloss_Beam,Rloss_Iso,Rloss_Albedo] = pvl_iam_martinruiz_components(SurfTilt,AOI,)|
%   |[Rloss_Beam,Rloss_Iso,Rloss_Albedo] = pvl_iam_martinruiz_components(SurfTilt,AOI,Rloss_Para)|
%
% Description
% This model calculates the reflection losses at the air/glass interface of
% a solar module separately for direct (beam), isotropic sky diffuse and 
% ground-reflected (albedo) irradiance. The model uses empirical equations 
% developed in [1].
%
% Inputs:
%   |SurfTilt| - a scalar or vector of surface tilt angles in decimal degrees.
%     If |SurfTilt| is a vector it must be of the same size as all other vector
%     inputs. |SurfTilt| must be >=0 and <=180. The tilt angle is defined as
%     degrees from horizontal (e.g. surface facing up = 0, surface facing
%     horizon = 90).
%   |AOI| - a scalar or vector of the angle of incidence of beam irradiance
%     in decimal degrees. If AOI is a vector it must be of the same size as all other vector
%     inputs. |AOI| must be >=0 and <=90. The angle of incidence can be calculated by
%     |pvl_getaoi|.
%   |Rloss_Para| - a three-element vector represents the parameters 
%     (in order, ar, c1, and c2) in the reflection models in Ref. [4].
%     By default a parameter set for a glass-faced silicon solar module, 
%     [ar = 0.16, cl = 4/(3*pi), c2 = -0.074], will be used.  
%
% Output:
%   |Rloss_Beam| - A column vector with the same number of elements as any input
%     vectors, which contains the reflection loss of the direct beam light. 
%   |Rloss_Iso| - A single value, the reflection loss of the isotropic diffuse irradiance. 
%   |Rloss_Albedo| - A singel value, the reflection loss of the ground-reflected light. 
%
% References
%   [1] Martín, N., Ruiz, J. M. 2005. Annual angular reflection losses in
%   PV modules. Progress in Photovoltaics: Research and Applications, 13(1), 75–84.

p = inputParser;
p.addRequired('SurfTilt', @(x) all(isnumeric(x) & x<=180 & x>=0 & isvector(x)));
p.addRequired('AOI', @(x) all(isnumeric(x) & x<=90 & x>=0 & isvector(x)));
p.addOptional('Rloss_Para', [0.16, 4/(3*pi), -0.074], @(x) isvector(x) && length(x)==3 );
p.parse(SurfTilt, AOI, varargin{:})

SurfTilt = p.Results.SurfTilt(:);
AOI = p.Results.AOI(:);
Rloss_Para= p.Results.Rloss_Para(:);

%Examin input size
VectorSizes = [numel(SurfTilt), numel(AOI)];
MaxVectorSize = max(VectorSizes);
if not(all((VectorSizes==MaxVectorSize) | (VectorSizes==1)))
    error(['Input parameters SurfTilt and AOI'...
        ' must either be scalars or vectors of the same length.']);
end

%Load the parameter
ar = Rloss_Para(1); 
c1 = Rloss_Para(2);
c2 = Rloss_Para(3);

%Equation 3a in [1]
Rloss_Beam = 1 - (1-exp(-cosd(AOI)/ar))./(1-exp(-1/ar));

%Equation 3c in [1]
Rloss_Iso = exp( -1/ar*(c1*(sind(SurfTilt)+(pi-SurfTilt/180*pi-sind(SurfTilt))./(1+cosd(SurfTilt))) + c2* (sind(SurfTilt)+(pi-SurfTilt/180*pi-sind(SurfTilt))./(1+cosd(SurfTilt))).^2));
%Equation 3b in [1]
Rloss_Albedo = exp( -1/ar*(c1*(sind(SurfTilt)+(SurfTilt/180*pi-sind(SurfTilt))./(1+cosd(SurfTilt))) + c2* (sind(SurfTilt)+(SurfTilt/180*pi-sind(SurfTilt))./(1+cosd(SurfTilt))).^2)); 

