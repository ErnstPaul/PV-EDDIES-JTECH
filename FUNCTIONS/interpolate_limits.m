function [in_out_vel]=interpolate_limits(pt_i,pt_j,km_i,km_j,dkm,vel,dir)
% interpolate_limits.m
%  interpolate_limits(pt_i,pt_j,km_i,km_j,dkm,vel,dir) interpolate velocity
%  magnitudes across the four extremes of a closed conotur.
%
%  - pt_i and pt_j are the location in kilometric coordinate of the curve
%      extreme under consideration;
%  - km_i and km_j are kilometric coordinate matrices (coordinate origin is 
%          the south-east corner of the area);
%  - dkm is the distance to and from the extreme where velocity magnitude
%          is interpolated;
%  - vel is the velocity magnitude field, used to check the increase in
%          velocity across the closed contour.
%  - dir is a string that identify the type of extreme (North 'N', East
%          'E', South 'S' or West 'W');
%
%  OUTPUT:
%  - in_ou_vel contains the velocity magnitude 'dkm' kilometers before and
%          after the curve extreme; the two are compared in 'max_curve.m'.
%
%-------------------------
%   Ver. 2.1 Oct.2012
%   Ver. 2.0 Jan.2012
%   Ver. 1.3 Apr.2011
%   Ver. 1.2 May.2010
%   Ver. 1.1 Dec.2009
%   Authors: Francesco Nencioli, francesco.nencioli@univ-amu.fr
%            Charles Dong, cdong@atmos.ucla.edu
%-------------------------
%
% Copyright (C) 2009-2012 Francesco Nencioli and Charles Dong
%
% This file is part of the Vector-Geometry Eddy Detection Algorithm.
%
% The Vector-Geometry Eddy Detection Algorithm is free software: 
% you can redistribute it and/or modify it under the terms of the 
% GNU General Public License as published by the Free Software Foundation, 
% either version 3 of the License, or (at your option) any later version.
% 
% The Vector-Geometry Eddy Detection Algorithm is distributed in the 
% hope that it will be useful, but WITHOUT ANY WARRANTY; without even 
% the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
% PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with the Vector-Geometry Eddy Detection Algorithm.  
% If not, see <http://www.gnu.org/licenses/>.
%
%=========================

%-------------------------------------------------------------------
% interpolate vel across the curve along different directions depending on 
% different extremes
switch (dir)
    case 'N'
        pts=[pt_i,pt_j-dkm;pt_i,pt_j];
    case 'S'
        pts=[pt_i,pt_j+dkm;pt_i,pt_j];
    case 'W'
        pts=[pt_i+dkm,pt_j;pt_i,pt_j];
    case 'E'
        pts=[pt_i-dkm,pt_j;pt_i,pt_j];
end
% interpolate vel
for j=1:length(pts(:,1))
	% Find i and j index of closest grid point to pts point
	aa=find(km_i<=pts(j,1) & km_j<=pts(j,2));
	minD=sqrt((km_i(aa)-pts(j,1)).^2+(km_j(aa)-pts(j,2)).^2);
	[J I]=ind2sub(size(km_i),aa(minD==min(minD)));
	% Closest point position in fraction of index
	if I==size(km_i,2) % To fix case where contour is out of domain 
		% (Matlab bug?)
		delta_i=1;
		I=I-1;
	else
		delta_i=(pts(j,1)-km_i(J,I))/(km_i(J,I+1)-km_i(J,I));
	end
	if J==size(km_j,1) % To fix case where contour is out of domain 
		% (Matlab bug?)
		delta_j=1;
		J=J-1;
	else
		delta_j=(pts(j,2)-km_j(J,I))/(km_j(J+1,I)-km_j(J,I));
	end
	% Velocity of pts point
	%         % (bi-linear interpolator)
	in_out_vel(j,1)=vel(J,I)*(1-delta_i)*(1-delta_j)+ ...
		vel(J,I+1)*delta_i*(1-delta_j)+ ...
		vel(J+1,I)*(1-delta_i)*delta_j+ ...
		vel(J+1,I+1)*delta_i*delta_j;
end
