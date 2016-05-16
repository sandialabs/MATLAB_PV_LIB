%% pvl_getaoi 
% Determine solar angle of incidence on a tilted surface.
%
%% Syntax
% |AOI = pvl_getaoi(SurfTilt, SurfAz, SunZen, SunAz)|
%   
%% Description
% Determines the angle of incidence in degrees between a surface normal and a
% vector pointed at the sun. The surface is defined by its tilt angle from
% horizontal and its azimuth pointing angle. The sun position is defined
% by the apparent (refraction corrected) sun zenith angle and the sun 
% azimuth angle.
%
%% Inputs
%%
% * *|SurfTilt|* - a scalar or vector of surface tilt angles in decimal degrees.
%     If |SurfTilt| is a vector it must be of the same size as all other vector
%     inputs. |SurfTilt| must be >=0 and <=180. The tilt angle is defined as
%     degrees from horizontal (e.g. surface facing up = 0, surface facing
%     horizon = 90)
% * *|SurfAz|* - a scalar or vector of surface azimuth angles in decimal degrees.
%     If |SurfAz| is a vector it must be of the same size as all other vector
%     inputs. |SurfAz| must be >=0 and <=360. The azimuth convention is defined
%     as degrees east of north (e.g. North = 0, East = 90, West = 270).
% * *|SunZen|* - a scalar or vector of apparent (refraction-corrected) zenith
%     angles in decimal degrees. If |SunZen| is a vector it must be of the
%     same size as all other vector inputs. |SunZen| must be >=0 and <=180.
% * *|SunAz|* - a scalar or vector of sun azimuth angles in decimal degrees.
%     If |SunAz| is a vector it must be of the same size as all other vector
%     inputs. |SunAz| must be >=0 and <=360. The azimuth convention is defined
%     as degrees east of north (e.g. North = 0, East = 90, West = 270).
%
%% Outputs
%%
% * *|AOI|* - A column vector with the same number of elements as any input
%     vectors, which contains the angle, in decimal degrees, between the 
%     surface normal vector and the sun beam vector. 
%
%% Example 1
% Sun at zenith
SurfTilt = 30;
SurfAz = 180;
SunZen = 0;
SunAz = 180;
AOI = pvl_getaoi(SurfTilt, SurfAz, SunZen, SunAz)
%% Example 2
% Sun low in southwestern sky 
SurfTilt = 30;
SurfAz = 180;
SunZen = 30;
SunAz = 225;
AOI = pvl_getaoi(SurfTilt, SurfAz, SunZen, SunAz)
%% References
% [1] D.L. King, J.A. Kratochvil, W.E. Boyson , 1997. Spectral and
%   Angle-of-Incidence Effects on Photovoltaic Modules and Solar Irradiance
%   Sensors. 26th IEEE Photovoltaic Specialists Conference.
%
%% See Also  
% <pvl_ephemeris_help.html |pvl_ephemeris|>
%%
% Copyright 2014 Sandia National Laboratories
