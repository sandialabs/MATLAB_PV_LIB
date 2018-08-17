function [TrkrTheta, AOI, SurfTilt, SurfAz] = pvl_singleaxis(SunZen, SunAz, Latitude, AxisTilt, AxisAzimuth, MaxAngle, varargin)
% PVL_SINGLEAXIS Determine the rotation angle of a 1 axis tracker, and sun incident angle on the tracked surface 
% 
% Syntax 
%   [TrkrTheta, AOI, SurfTilt, SurfAz]= pvl_singleaxis(SunZen, SunAz, Latitude, AxisTilt, AxisAzimuth, MaxAngle)
%   [TrkrTheta, AOI, SurfTilt, SurfAz]= pvl_singleaxis(SunZen, SunAz, Latitude, AxisTilt, AxisAzimuth, MaxAngle, Backtrack)
%   [TrkrTheta, AOI, SurfTilt, SurfAz]= pvl_singleaxis(SunZen, SunAz, Latitude, AxisTilt, AxisAzimuth, MaxAngle, Backtrack, GCR)
% 
% Description
%   Determine the rotation angle of a single axis tracker using the
%   equations in [1] when given a particular sun zenith and azimuth angle.
%   Backtracking may be specified, and if so, a ground coverage ratio is 
%   required.
%
%   Rotation angle is determined in a panel-oriented coordinate system.
%   The tracker azimuth AxisAzimuth defines the positive y-axis; the positive
%   x-axis is 90 degress clockwise from the y-axis and parallel to the
%   earth surface, and the positive z-axis is normal and oriented towards the sun.
%   Rotation angle TrkrTheta indicates tracker position relative to horizontal:
%   TrkrTheta = 0 is horizontal, and positive TrkrTheta is a clockwise rotation
%   around the y axis in the x, y, z coordinate system.  For example, if tracker azimuth 
%   AxisAzimuth is 180 (oriented south), TrkrTheta = 30 is a rotation of 
%   30 degrees towards the west, and TrkrTheta = -90 is a rotation to the 
%   vertical plane facing east.
% 
% Inputs:
%   SunZen - a scalar or vector of apparent (refraction-corrected) zenith
%     angles in decimal degrees. If SunZen is a vector it must be of the
%     same size as all other vector inputs. SunZen must be >=0 and <=180.
%   SunAz - a scalar or vector of sun azimuth angles in decimal degrees.
%     If SunAz is a vector it must be of the same size as all other vector
%     inputs. SunAz must be >=0 and <=360. The Azimuth convention is defined
%     as degrees East of North (e.g. North = 0, East = 90, West = 270).
%   Latitude - a scalar value denoting which hempisphere the tracker is
%     in. The exact latitude is NOT required, any positive number denotes
%     the northern hemisphere, any negative number denotes the southern
%     hemisphere, a value of 0 is assumed to be northern hemisphere.
%   AxisTilt - a scalar value denoting the tilt of the axis of rotation
%     (i.e, the y-axis defined by AxisAzimuth) with respect to horizontal, 
%     in decimal degrees. AxisTilt must be >=0 and <=180.
%   AxisAzimuth - a scalar value denoting the compass direction along which
%     the axis of rotation lies, in decimal degrees. Again, the convention 
%     is defined as degrees East of North (e.g. North = 0, East = 90, 
%     West = 270). AxisAzimuth must be >=0 and <=360.
%   MaxAngle - a scalar value denoting the maximum rotation angle, in
%     decimal degrees, of the one-axis tracker from its horizontal position
%     (horizontal if AxisTilt = 0). MaxAngle must be <=180 and >=0. A
%     MaxAngle of 90 degrees allows the tracker to rotate to a vertical
%     position to point the panel towards a horizon.  
%     MaxAngle of 180 degrees allows for full rotation.
%   Backtrack - a scalar value denoting whether the tracker has the
%     capability to "backtrack" to avoid row-to-row shading. A value of 0
%     denotes no backtrack capability. Any other value denotes backtrack
%     capability. If Backtrack is not provided, the default value is 0 (no
%     backtrack).
%   GCR - a scalar value denoting the ground coverage ratio of a tracker
%     system which utilizes backtracking; i.e. the ratio between the PV
%     array surface area to total ground area. A tracker system with modules 2
%     meters wide, centered on the tracking axis, with 6 meters between the
%     tracking axes has a GCR of 2/6=0.333. If GCR is not provided, a GCR
%     of 2/7 is default. GCR must be <=1.
% 
% Output:
%   TrkrTheta - The rotation angle (Theta) of the tracker.  
%     TrkrTheta = 0 is horizontal, and positive rotation angles are
%     clockwise.
%   AOI - The angle-of-incidence of direct irradiance onto the
%     rotated panel surface.
%   SurfTilt - The angle between the panel surface and the earth
%     surface, accounting for panel rotation.
%   SurfAz - The azimuth of the rotated panel, determined by 
%     projecting the vector normal to the panel's surface to the earth's
%     surface.
%
% References:
%   [1] Lorenzo, E et al., 2011, "Tracking and back-tracking", Prog. in 
%   Photovoltaics: Research and Applications, v. 19, pp. 747-753.
%
% See Also: PVL_EPHEMERIS PVL_SPA
%
%
p = inputParser;
p.addRequired('SunZen', @(x) all(isnumeric(x) & x>=0 & x<=180) & isvector(x));
p.addRequired('SunAz', @(x) all(isnumeric(x) & x>=0 & x<=360) & isvector(x));
p.addRequired('Latitude', @(x) (isnumeric(x) & isscalar(x)));
p.addRequired('AxisTilt', @(x) (isnumeric(x) & isscalar(x) & x>=0 & x<=90));
p.addRequired('AxisAzimuth', @(x) (isnumeric(x) & isscalar(x) & x>=0 & x<=360));
p.addRequired('MaxAngle', @(x) (isnumeric(x) & isscalar(x) & x<=180 & x>=0));
p.addOptional('BackTrack', 0, @(x) (isscalar(x)));
p.addOptional('GCR',1/3.5, @(x) (isnumeric(x) & isscalar(x) & x<=1));
p.parse(SunZen, SunAz, Latitude, AxisTilt, AxisAzimuth, MaxAngle, varargin{:})

if Latitude<0
    msgbox('This code doesn''t work for the southern hemisphere');
    return
end;

BackTrack = p.Results.BackTrack;
GCR = p.Results.GCR;
SunZen = SunZen(:);
SunAz = SunAz(:);
% Calculate sun position x, y, z using coordinate system as in [1], Eq 2.
% Positive y axis is oriented parallel to earth surface along tracking axis 
% (for the purpose of illustration, assume y is oriented to the south);
% positive x axis is orthogonal, 90 deg clockwise from y-axis, and parallel
% to the earth's surface (if y axis is south, x axis is west); 
% positive z axis is normal to x,y axes, pointed upward.
% Equations in [1] assume solar azimuth is relative to reference vector
% pointed south, with clockwise positive.  Here, the input solar azimuth 
% is degrees East of North, i.e., relative to a reference vector pointed 
% north with clockwise positive.  Rotate sun azimuth to coordinate system as in [1]
% to calculate sun position.
Az = SunAz - 180;
El = 90-SunZen;
x = cosd(El).*sind(Az);
y = cosd(El).*cosd(Az);
z = sind(El);

% translate array azimuth from compass bearing to [1] coord system
AxisAz = AxisAzimuth - 180;

% translate input array tilt angle axistilt to [1] coordinate system.  In 
% [1] coordinates, axistilt is a rotation about the x-axis.  For a system 
% with array azimuth (y-axis) oriented south, the x-axis is oriented west,
% and a positive axistilt is a counterclockwise rotation, i.e, lifting the 
% north edge of the panel.  Thus, in [1] coordinate system, in the northern 
% hemisphere a positive axistilt indicates a rotation toward the equator, 
% whereas in the southern hemisphere rotation toward the equator is 
% indicated by axistilt<0.  Here, the input axistilt is always positive and
% is a rotation toward the equator.

% Calculate sun position (xp, yp, zp) in panel-oriented coordinate system: 
% positive y-axis is oriented along tracking axis at panel tilt;
% positive x-axis is orthogonal, clockwise, parallel to earth surface;
% positive z-axis is normal to x-y axes, pointed upward.  
% Calculate sun position (xp,yp,zp) in panel coordinates using [1] Eq 11
xp = x.*cosd(AxisAz) - y.*sind(AxisAz);
yp = x.*cosd(AxisTilt).*sind(AxisAz) + y.*cosd(AxisTilt).*cosd(AxisAz) - z.*sind(AxisTilt);
% note that equation for yp (y' in Eq. 11 of Lorenzo et al 2011) is
% corrected, after conversation with paper's authors
zp = x.*sind(AxisTilt).*sind(AxisAz) + y.*sind(AxisTilt).*cosd(AxisAz) + z.*cosd(AxisTilt);

% The ideal tracking angle wid is the rotation to place the sun position 
% vector (xp, yp, zp) in the (y, z) plane; i.e., normal to the panel and 
% containing the axis of rotation.  wid = 0 indicates that the panel is 
% horizontal.  Here, our convention is that a clockwise rotation is 
% positive, to view rotation angles in the same frame of reference as 
% azimuth.  For example, for a system with tracking axis oriented south, 
% a rotation toward the east is negative, and a rotation to the west is 
% positive.

% filter to avoid undefined inverse tangent
tmp(xp~=0) = atand(zp./xp);  % angle from x-y plane to projection of sun vector onto x-z plane
tmp(xp==0 & zp>=0) = 90;    % fill in when atan is undefined
tmp(xp==0 & zp<0) = -90;    % fill in when atan is undefined
tmp=tmp(:);                  % ensure tmp is a column vector
% Obtain wid by translating tmp to convention for rotation angles.
% Have to account for which quadrant of the x-z plane in which the sun 
% vector lies.  Complete solution here but probably not necessary to 
% consider QIII and QIV.
wid(xp>=0 & zp>=0) =  90 - tmp(xp>=0 & zp>=0);  % QI
wid(xp<0  & zp>=0) = -90 - tmp(xp<0  & zp>=0);  % QII
wid(xp<0  & zp<0)  = -90 - tmp(xp<0  & zp<0);   % QIII
wid(xp>=0 & zp<0)  =  90 - tmp(xp>=0 & zp<0);   % QIV
wid=wid(:);                  % ensure wid is a column vector

% filter for sun above panel horizon)
u = zp>0;

% apply limits to ideal rotation angle
wid(~u) = 0;  % set horizontal if zenith<0, sun is below panel horizon

% Account for backtracking; modified from [1] to account for rotation
% angle convention being used here.
if (BackTrack ~= 0)
    Lew = 1/GCR;
    temp = min(Lew.*cosd(wid),1);
    wc = acosd(temp);   % backtrack angle; always positive (acosd returns values between 0 and 180)
    v = wid<0;
    widc(~v) = wid(~v) - wc(~v); % Eq 4 applied when wid in QI
    widc(v) = wid(v) + wc(v);    % Eq 4 applied when wid in QIV
else
    widc = wid;
end

TrkrTheta(u) = widc(u);
TrkrTheta(~u) = 0;    % set to zero when sun is below panel horizon
TrkrTheta = TrkrTheta(:);   % ensure column vector format
TrkrTheta(TrkrTheta > MaxAngle) = MaxAngle;
TrkrTheta(TrkrTheta < -MaxAngle) = -MaxAngle;

% calculate normal vector to panel in panel-oriented x, y, z coordinates
% y-axis is axis of tracker rotation.  TrkrTheta is a compass angle
% (clockwise is positive) rather than a trigonometric angle.

Norm = [sind(TrkrTheta) zeros(length(TrkrTheta),1) cosd(TrkrTheta)];

% sun position in vector format in panel-oriented x, y, z coordinates
P = [xp yp zp];

% calculate angle-of-incidence on panel
AOI = acosd(abs(dot(Norm,P,2)));
AOI(~u) = 0;    % set to zero when sun is below panel horizon

% calculate panel elevation SurfEl and azimuth SurfAz in a coordinate system where the
% panel elevation is the angle from horizontal, and the panel azimuth is
% the compass angle (clockwise from north) to the projection of the panel's
% normal to the earth's surface.  These outputs are provided for
% convenience and comparison with other PV software which use these angle
% conventions.

% project normal vector to earth surface.  First rotate 
% about x-axis by angle -AxisTilt so that y-axis is also parallel to earth 
% surface, then project.
Rot_x = [1 0 0; 0 cosd(-AxisTilt) -sind(-AxisTilt); 0 sind(-AxisTilt) cosd(-AxisTilt)]; %rotation matrix
temp = Rot_x * Norm'; % temp contains the normal vector expressed in earth-surface coordinates (z normal to surface, y aligned with tracker axis parallel to earth)
temp = temp';  % ensure column format
projNorm = [temp(:,1) temp(:,2) zeros(size(temp(:,3)))]; % projection to plane tangent to earth surface, in earth surface coordinates
tempnorm = sqrt(sum(temp.^2,2));
projNormnorm = sqrt(sum(projNorm.^2,2));

SurfAz = 0.*TrkrTheta;
% calculation of SurfAz
SurfAz = atand(projNorm(:,2)./projNorm(:,1));

% clean up atan when x-coord is zero
SurfAz(projNorm(:,1)==0 & projNorm(:,2)>0) =  90;
SurfAz(projNorm(:,1)==0 & projNorm(:,2)<0) =  -90;
% clean up atan when y-coord is zero
SurfAz(projNorm(:,2)==0 & projNorm(:,1)>0) =  0;
SurfAz(projNorm(:,2)==0 & projNorm(:,1)<0) = 180;
% correct for QII and QIII
SurfAz(projNorm(:,1)<0 & projNorm(:,2)>0) =  SurfAz(projNorm(:,1)<0 & projNorm(:,2)>0) + 180;  % QII
SurfAz(projNorm(:,1)<0 & projNorm(:,2)<0) =  SurfAz(projNorm(:,1)<0 & projNorm(:,2)<0) + 180;  % QIII

% at this point SurfAz contains angles between -90 and +270, where 0 is
% along the positive x-axis, the y-axis is in the direction of the tracker
% azimuth, and positive angles are rotations from the positive x axis towards
% the positive y-axis.
% Adjust to compass angles (clockwise rotation from 0 along the positive y-axis)
SurfAz(SurfAz<=90) = 90 - SurfAz(SurfAz<=90);
SurfAz(SurfAz>90) = 450 - SurfAz(SurfAz>90);
% finally rotate to align y-axis with true north
if Latitude>0
    SurfAz = SurfAz - AxisAzimuth;
else
    SurfAz = SurfAz - AxisAzimuth - 180;
end;
SurfAz(SurfAz<0) = 360+ SurfAz(SurfAz<0);

divisor = (round(tempnorm.*projNormnorm.*10000))./10000;
dividend = (round(dot(temp,projNorm,2).*10000))./10000;


SurfTilt = 90 - acosd(dividend./divisor);
SurfTilt = SurfTilt(:);

end

