function [Front_Irradiance,Rear_Irradiance] = pvl_Purdue_bifacial_irradiance(SurfTilt, SurfAz, EtoH, Albedo, DHI, DNI, HExtra, SunZen, SunAz, AM, varargin)

% |pvl_Purdue_bifacial_irradiance| calculates the irradiance on the
% front and rear sides of a bifacial solar module while fully accounting for
% the self-shading losses.
%
% Syntax
%   |pvl_Purdue_bifacial_irradiance(SurfTilt, SurfAz, EtoH, DHI, DNI, Albedo, HExtra, SunZen, SunAz, AM)|
%   |pvl_Purdue_bifacial_irradiance(SurfTilt, SurfAz, EtoH, DHI, DNI, Albedo, HExtra, SunZen, SunAz, AM, model)|
%   |pvl_Purdue_bifacial_irradiance(SurfTilt, SurfAz, EtoH, DHI, DNI, Albedo, HExtra, SunZen, SunAz, AM, model, Rloss)|
%   |pvl_Purdue_bifacial_irradiance(SurfTilt, SurfAz, EtoH, DHI, DNI, Albedo, HExtra, SunZen, SunAz, AM, model, Rloss, Rloss_Para)|
%
% Description
%   The Purdue Bifacial irradiance model [1] simulates the total irradiance including 
%   direct, diffuse, and albedo light, on both the front and rear sides of
%   a bifacial solar module. This model applies an analytical view-factor 
%   based approach to explicitly account for the self-shading losses of 
%   albedo light due to (1) direct blocking of direct and circumsolar diffuse
%   light and 2) sky masking of isotropic diffuse light onto the ground. 
%   This model also incorporates an optional reflection-loss model [4]. 
%   This model has been validated against data spanning from Africa, Europe,
%   Asia, and North America (please refer [1] for more detail).
%
% Inputs:
%   |SurfTilt| - a scalar or vector of surface tilt angles in decimal degrees.
%     The tilt angle is defined as degrees from horizontal (e.g. surface facing
%     up = 0, surface facing horizon = 90). |SurfTilt| must be >=0 and <=180.
%     If |SurfTilt| is a vector it must be of the same size as all other vector inputs.  
%   |SurfAz| - a scalar or vector of surface azimuth angles in decimal degrees.
%     If |SurfAz| is a vector it must be of the same size as all other vector
%     inputs. |SurfAz| must be >=0 and <=360. The Azimuth convention is defined
%     as degrees east of north (e.g. North = 0, East = 90, West = 270).
%   |EtoH| - a scalar or vector of the ratio of module elevation(E) to module height(H).
%     Module height is the module dimension not parallel to the ground.
%     If |EtoH| is a vector it must be of the same size as all other vector
%     inputs. |EtoH| must be >=0.
%   |Albedo| - a scalar or vector of groud albedo coefficient.
%     If |Albedo| is a vector it must be of the same size as all other vector
%     inputs. |Albedo| must be >=0 and <=1.
%   |DHI| - a scalar or vector of diffuse horizontal irradiance in W/m^2. 
%     If |DHI| is a vector it must be of the same size as all other vector inputs.
%     |DHI| must be >=0.
%   |DNI| - a scalar or vector of direct normal irradiance in W/m^2. If
%   |DNI| is a vector it must be of the same size as all other vector inputs.
%     |DNI| must be >=0.
%   |HExtra| - a scalar or vector of extraterrestrial normal irradiance in
%     W/m^2. If |HExtra| is a vector it must be of the same size as
%     all other vector inputs. |HExtra| must be >=0.
%   |SunZen| - a scalar or vector of apparent (refraction-corrected) zenith
%     angles in decimal degrees. If |SunZen| is a vector it must be of the
%     same size as all other vector inputs. |SunZen| must be >=0 and <=180.
%   |SunAz| - a scalar or vector of sun azimuth angles in decimal degrees.
%     If |SunAz| is a vector it must be of the same size as all other vector
%     inputs. |SunAz| must be >=0 and <=360. The Azimuth convention is defined
%     as degrees east of north (e.g. North = 0, East = 90, West = 270).
%   |AM| - a scalar or vector of relative (not pressure-corrected) airmass 
%     values. If |AM| is a vector it must be of the same size as all other 
%     vector inputs. |AM| must be >=0.
%   |model| - a character string which selects the desired set of Perez
%     coefficients. If model is not provided as an input, the default,
%     '1990' will be used.
%     All possible model selections are: 
%       '1990', 'allsitescomposite1990' (same as '1990'),
%       'allsitescomposite1988', 'sandiacomposite1988',
%       'usacomposite1988', 'france1988', 'phoenix1988',
%       'elmonte1988', 'osage1988', 'albuquerque1988',
%       'capecanaveral1988', or 'albany1988'
%   |Rloss| - a character string which determines the inclusion of reflection 
%     loss model. By default, 'Yes' will be used. If 'No' is input, reflection 
%     loss will be neglected.
%   |Rloss_Para| - a three-element vector represents the parameters 
%     (in order, ar, c1, and c2) in the reflection models in Ref. [4].
%     By default a parameter set for a glass-faced silicon solar module, 
%     [ar = 0.16, cl = 4/(3*pi), c2 = -0.074], will be used.  
%
% Output:
%   |Front_Irradiance| - the total irradiance including direct beam, diffuse,
%     and albedo light on the front side. |Front_Irradiance| is a column vector
%     of the same size as the input vector(s).
%   |Rear_Irradiance| - the total irradiance includig direct beam, diffuse,
%     and albedo light on the rear side.  |Rear_Irradiance| is a column vector
%     of the same size as the input vector(s).

% References
%   [1] Sun, X., Khan, M. R., Alam, M. A., 2018. Optimization and performance 
%   of bifacial solar modules: A global perspective. Applied Energy 212, pp. 1601-1610.
%   [2] Khan, M. R., Hanna, A., Sun, X., Alam, M. A., 2017. Vertical bifacial solar farms:
%   Physics, design, and global optimization. Applied Energy, 206, 240–248.
%   [3] Duffie, J. A., Beckman, W. A. 2013. Solar Engineering of Thermal Processes (4th Editio). 
%   Wiley.
%   [4] Martín, N., Ruiz, J. M. 2005. Annual angular reflection losses in
%   PV modules. Progress in Photovoltaics: Research and Applications, 13(1), 75–84.
%
% See also |pvl_perez| |pvl_Purdue_albedo_model|
% |pvl_iam_martinruiz_components|
%
% Notes: pvl_Purdue_bifacial_irradiance contributed by Xingshu Sun of Purdue
% University, 2018.

%% Process Inputs
%parse parameters
p=inputParser;
p.addRequired('SurfTilt', @(x) (isnumeric(x) && all(x<=180) && all(x>=0) && isvector(x)));
p.addRequired('SurfAz', @(x) isnumeric(x) && all(x<=360) && all(x>=0) && isvector(x));
p.addRequired('EtoH', @(x) isnumeric(x) && all(x>=0) && isvector(x));
p.addRequired('Albedo', @(x) isnumeric(x) && all(x<=1) && all(x>=0) && isvector(x));
p.addRequired('DHI', @(x) (isnumeric(x) && isvector(x) && all((x>=0) | isnan(x))));
p.addRequired('DNI', @(x) isnumeric(x) && isvector(x) && all((x>=0) | isnan(x)));
p.addRequired('HExtra', @(x) isnumeric(x) && isvector(x) && all((x>=0) | isnan(x)));
p.addRequired('SunZen', @(x) isnumeric(x) && all(x<=180) && all((x>=0) | isnan(x)) && isvector(x));
p.addRequired('SunAz', @(x) (isnumeric(x) && all(x<=360) && all((x>=0) | isnan(x)) && isvector(x)));
p.addRequired('AM', @(x) (all(((isnumeric(x) & x>=0) | isnan(x))) & isvector(x)));
p.addOptional('model', '1990', @(x) ischar(x));
p.addOptional('Rloss', 'Yes', @(x) ischar(x));
p.addOptional('Rloss_Para', [0.16, 4/(3*pi), -0.074], @(x) isvector(x) && length(x)==3 );
p.parse(SurfTilt, SurfAz, EtoH, Albedo, DHI, DNI, HExtra, SunZen, SunAz, AM, varargin{:});


SurfTilt = p.Results.SurfTilt(:);
SurfAz = p.Results.SurfAz(:);
EtoH = p.Results.EtoH(:);
Albedo = p.Results.Albedo(:);
DHI = p.Results.DHI(:);
DNI = p.Results.DNI(:);
HExtra = p.Results.HExtra(:);
SunZen = p.Results.SunZen(:);
SunAz = p.Results.SunAz(:);
AM = p.Results.AM(:);
model = p.Results.model;
Rloss = p.Results.Rloss;
Rloss_Para = p.Results.Rloss_Para;


%Examin input size
VectorSizes = [numel(SurfTilt), numel(SurfAz), numel(DHI), numel(DNI), ...
    numel(HExtra), numel(SunZen), numel(SunAz), numel(AM)];
MaxVectorSize = max(VectorSizes);
if not(all((VectorSizes==MaxVectorSize) | (VectorSizes==1)))
    error(['Input parameters SurfTilt, SurfAz, EtoH, DHI, DNI, Albedo, HExtra, SunZen, SunAz, AM'...
        ' must either be scalars or vectors of the same length.']);
end

%% Irradiance calculation
%%%%%%%%%%%%%%%%%%%%Front Side%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                               
%Direct beam 
AOI_Front = pvl_getaoi(SurfTilt, SurfAz, SunZen, SunAz);
AOI_Front((AOI_Front>90)|(AOI_Front<0))=90;
IB_Front = DNI .* cosd(AOI_Front); 
%Account for reflection loss
if strcmp(Rloss,'Yes')
    [Rloss_Beam_Front, Rloss_Iso_Front, Rloss_Albedo_Front] = pvl_iam_martinruiz_components(SurfTilt,AOI_Front,Rloss_Para);
    % for horizon brightening
    AOI_Hor_Front = pvl_getaoi(SurfTilt, SurfAz, 90, SurfAz);
    AOI_Hor_Front((AOI_Hor_Front>90)|(AOI_Hor_Front<0))=90;
    Rloss_Hor_Front = pvl_iam_martinruiz_components(SurfTilt,AOI_Hor_Front,Rloss_Para);
else
    Rloss_Beam_Front = 0;
    Rloss_Iso_Front = 0;
    Rloss_Albedo_Front = 0;
    Rloss_Hor_Front = 0;
end

IB_Front  = IB_Front .* (1-Rloss_Beam_Front);
IB_Front(IB_Front<0) = 0;
                                         
%Sky diffuse
[~,ID_Iso_Front,ID_Cir_Front,ID_Hor_Front] = pvl_perez(SurfTilt, SurfAz,...
    DHI, DNI,HExtra, SunZen, SunAz, AM ,model); %Perez Diffuse  

ID_Iso_Front = ID_Iso_Front .* (1-Rloss_Iso_Front); 
ID_Cir_Front = ID_Cir_Front .* (1-Rloss_Beam_Front);
ID_Hor_Front = ID_Hor_Front .* (1-Rloss_Hor_Front);
    
%Albedo light
I_Alb_Front = pvl_Purdue_albedo_model(SurfTilt, SurfAz, EtoH, Albedo, DHI, DNI, HExtra, SunZen, SunAz, AM ,model);
I_Alb_Front = I_Alb_Front .* (1-Rloss_Albedo_Front);
I_Alb_Front(I_Alb_Front<0) =0;

%Sum up the front-side irradiance
Front_Irradiance = IB_Front + I_Alb_Front + ID_Iso_Front+ ID_Cir_Front+ ID_Hor_Front;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%Rear Side%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Define angle for the rear side
SurfTilt_Rear = 180 - SurfTilt;
SurfAz_Rear = SurfAz+180;
if SurfAz_Rear >= 360
   SurfAz_Rear = SurfAz_Rear - 360;
end                                                          
                
%Direct beam
AOI_Rear = pvl_getaoi(SurfTilt_Rear, SurfAz_Rear, SunZen, SunAz);
AOI_Rear((AOI_Rear>90)|(AOI_Rear<0))=90;
IB_Rear = DNI .* cosd(AOI_Rear); 
%Account for reflection loss
if strcmp(Rloss,'Yes')
    [Rloss_Beam_Rear,Rloss_Iso_Rear,Rloss_Albedo_Rear] = pvl_iam_martinruiz_components(SurfTilt_Rear,AOI_Rear,Rloss_Para);
    %Horizon Brightening
    AOI_Hor_Rear = pvl_getaoi(SurfTilt_Rear, SurfAz_Rear, 90, SurfAz_Rear);
    AOI_Hor_Rear((AOI_Hor_Rear>90)|(AOI_Hor_Rear<0))=90;
    Rloss_Hor_Rear = pvl_iam_martinruiz_components(SurfTilt_Rear,AOI_Hor_Rear,Rloss_Para);
else
    Rloss_Beam_Rear = 0;
    Rloss_Iso_Rear = 0;
    Rloss_Albedo_Rear = 0;
    Rloss_Hor_Rear = 0;
end

IB_Rear = IB_Rear .* (1-Rloss_Beam_Rear); 
IB_Rear(IB_Rear<0) = 0;
                                               
%Sky diffuse light
[~,ID_Iso_Rear,ID_Cir_Rear,ID_Hor_Rear] = pvl_perez(SurfTilt_Rear, SurfAz_Rear,...
   DHI, DNI,HExtra, SunZen, SunAz, AM, model ); %Perez Diffuse

ID_Iso_Rear = ID_Iso_Rear .* (1-Rloss_Iso_Rear);        
ID_Cir_Rear = ID_Cir_Rear .* (1-Rloss_Beam_Rear);
ID_Hor_Rear = ID_Hor_Rear .* (1-Rloss_Hor_Rear);  

%Albedo light
I_Alb_Rear = pvl_Purdue_albedo_model(SurfTilt_Rear, SurfAz_Rear, EtoH, Albedo, DHI, DNI, HExtra, SunZen, SunAz, AM ,model);
I_Alb_Rear = I_Alb_Rear .* (1-Rloss_Albedo_Rear);
I_Alb_Rear(I_Alb_Rear<0) = 0;

%Sum up the rear-side irradiance
Rear_Irradiance =  IB_Rear + I_Alb_Rear + ID_Iso_Rear+ ID_Cir_Rear+ ID_Hor_Rear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
                  