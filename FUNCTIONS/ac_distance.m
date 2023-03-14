function [rng, az] = ac_distance(lat1,lon1,lat2,lon2)
%DISTANCE  Distance between points on sphere or ellipsoid
%
%   [ARCLEN, AZ] = DISTANCE(LAT1,LON1,LAT2,LON2) computes the lengths,
%   ARCLEN, of the great circle arcs connecting pairs of points on the
%   surface of a sphere. In each case, the shorter (minor) arc is assumed.
%   The function can also compute the azimuths, AZ, of the second point in
%   each pair with respect to the first (that is, the angle at which the
%   arc crosses the meridian containing the first point).  The input
%   latitudes and longitudes, LAT1, LON1, LAT2, LON2, can be scalars or
%   arrays of equal size and must be expressed in degrees. ARCLEN is
%   expressed in degrees of arc and will have the same size as the input
%   arrays.  AZ is measured clockwise from north, in units of degrees.
%   When given a combination of scalar and array inputs, the scalar inputs
%   are automatically expanded to match the size of the arrays.
%
%   [ARCLEN, AZ] = DISTANCE(LAT1,LON1,LAT2,LON2,ELLIPSOID) computes
%   geodesic arc length and azimuth assuming that the points lie on the
%   reference ellipsoid defined by the input ELLIPSOID.  ELLIPSOID is a
%   reference ellipsoid (oblate spheroid) object, a reference sphere
%   object, or a vector of the form [semimajor_axis, eccentricity].  The
%   output, ARCLEN, is expressed in the same length units as the semimajor
%   axis of the ellipsoid.
%
%   [ARCLEN, AZ] = DISTANCE(LAT1,LON1,LAT2,LON2,UNITS) uses the input
%   string UNITS to define the angle unit of the outputs ARCLEN and AZ and
%   the input latitude-longitude coordinates.  UNITS may equal 'degrees'
%   (the default value) or 'radians'.
%
%   [ARCLEN, AZ] = DISTANCE(LAT1,LON1,LAT2,LON2,ELLIPSOID,UNITS) uses the
%   UNITS string to specify the units of the latitude-longitude
%   coordinates, but the output range has the same units as the
%   semimajor axis of the ellipsoid.
%
%   [ARCLEN, AZ] = DISTANCE(TRACK,...) uses the input string TRACK to
%   specify either a great circle/geodesic or a rhumb line arc. If TRACK
%   equals 'gc' (the default value), then great circle distances are
%   computed on a sphere and geodesic distances are computed on an
%   ellipsoid. If TRACK equals 'rh', then rhumb line distances are computed
%   on either a sphere or ellipsoid.
%
%   [ARCLEN, AZ] = DISTANCE(PT1,PT2) accepts N-by-2 coordinate arrays
%   PT1 and PT2 such that PT1 = [LAT1 LON1] and PT2 = [LAT2 LON2] where
%   LAT1, LON1, LAT2, and LON2 are column vectors.  It is equivalent to
%   ARCLEN = DISTANCE(PT1(:,1),PT1(:,2),PT2(:,1),PT2(:,2)).
%
%   [ARCLEN, AZ] = DISTANCE(PT1,PT2,ELLIPSOID)
%   [ARCLEN, AZ] = DISTANCE(PT1,PT2,UNITS),
%   [ARCLEN, AZ] = DISTANCE(PT1,PT2,ELLIPSOID,UNITS) and
%   [ARCLEN, AZ] = DISTANCE(TRACK,PT1,...) are all valid calling forms.
%
%   Remark on Computing Azimuths
%   ----------------------------
%   Note that when both distance and azimuth are required for the same
%   point pair(s), it's more efficient to compute both with a single
%   call to DISTANCE.  That is, use:
%
%       [ARCLEN, AZ] = DISTANCE(...);
%
%   rather than its slower equivalent:
%
%       ARCLEN = DISTANCE(...);
%       AZ = AZIMUTH(...);
%
%   Remark on Output Units
%   ----------------------
%   To express the output ARCLEN as an arc length expressed in either
%   degrees or radians, omit the ELLIPSOID argument.  This is possible only
%   on a sphere.  If ELLIPSOID is supplied, ARCLEN is expressed in the same
%   length units as the semimajor axis of the ellipsoid.  Specify ELLIPSOID
%   as [R 0] to compute ARCLEN as a distance on a sphere of radius R, with
%   ARCLEN having the same units as R.
%
%   Remark on Eccentricity
%   ----------------------
%   Geodesic distances on an ellipsoid are valid only for small
%   eccentricities typical of the Earth (e.g., 0.08 or less).
%
%   Remark on Long Geodesics
%   ------------------------
%   Distance calculations for geodesics degrade slowly with increasing
%   distance and may break down for points that are nearly antipodal,
%   and/or when both points are very close to the Equator.  In addition,
%   for calculations on an ellipsoid, there is a small but finite input
%   space, consisting of pairs in which both the points are nearly
%   antipodal AND both points fall close to (but not precisely on) the
%   Equator. In this case, a warning is issued and both ARCLEN and AZ
%   are set to NaN for the "problem pairs."
%
%   See also AZIMUTH, RECKON.

% Copyright 1996-2011 The MathWorks, Inc.
% $Revision: 1.10.4.12 $  $Date: 2011/09/19 17:44:44 $


lat1 = (pi/180) * lat1;
lat2 = (pi/180) * lat2;
lon1 = (pi/180) * lon1;
lon2 = (pi/180) * lon2;

% Start with spherical approximation even if using an ellipsoid
a = sin((lat2-lat1)/2).^2 + cos(lat1) .* cos(lat2) .* sin((lon2-lon1)/2).^2;
rng = 2 * atan2(sqrt(a),sqrt(1 - a));
%rng = reshape(rng,insize);
rng = 6371 * rng;


if nargout >1
    % Azimuth on a sphere
    
    
    az = atan2(cos(lat2) .* sin(lon2-lon1),cos(lat1) .* sin(lat2) - sin(lat1) .* cos(lat2) .* cos(lon2-lon1));
    
    % Azimuths are undefined at the poles, so we choose a convention: zero at
    % the north pole and pi at the south pole.
    az(lat1 <= -pi/2) = 0;
    az(lat2 >=  pi/2) = 0;
    az(lat2 <= -pi/2) = pi;
    az(lat1 >=  pi/2) = pi;
    
    % Ensure azimuth in the range [0 2*pi).
    az = mod(az, 2*pi);
    az = az*180/pi;
end
