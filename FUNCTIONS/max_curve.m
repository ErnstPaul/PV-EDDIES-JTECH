function [eddy_lim large]=max_curve(km_i,km_j,field,km_ci,km_cj,vel)
% max_curve.m
%  max_curve(km_i,km_j,field,km_ci,km_cj,vel) computes the eddy shape
%  defined as the largest closed contour of the streamfunction (PSI) field,
%  across which velocity magnitude increases.
%  The increase of velocity is checked at the four extremes of the closed
%  contour (northermost, easternmost, southermost and westernmost points).
%  If velocity never increases across a closed contour, the shape is
%  defined by simply the largest closed contour around the center.
%
%  - km_i and km_j are kilometric coordinate matrices (coordinate origin is 
%          the south-east corner of the area);
%  - field is the PSI field, from where the eddy shape is derived;
%  - km_ci and km_cj are the kilometric coordinates of the detected eddy
%          center;
%  - vel is the velocity magnitude field, used to check the increase in
%          velocity across the closed contour.
%
%  OUTPUT:
%  - eddy_lim is a 2xn array containing the poisition in kilometric 
%          coordinates of the n vertices that define eddy shape;
%         (first row x, second y positions)
% - large is the flag which will be used in 'eddy_dim.m', that indicates if
%          the shape is normally defined (0), or if it's defined as simply
%          the largest closed contour around the center (1);
%
% 'interpolate limits' is used to derive velocity magnitude across the four
% extremes of a closed contour. Check the documentation in 
% 'interpolate limits.m' for further details.
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

% initialize the flag
large=0;
% define the distance from the contour extremes at which velocity magnitude 
% is interpolated
[c_j, c_i]=find(km_j==km_cj & km_i==km_ci);
dkm=0.05 * ...      % distance is set in km
	((km_i(c_j,c_i+1)-km_i(c_j,c_i-1))/2 + ...
	(km_j(c_j+1,c_i)-km_j(c_j-1,c_i))/2)/2;
% compute contourlines of the streamfunction field (100 contours)
C=contourc(km_i(1,:),km_j(:,1),field,100);
% C is a vector conatining all the the PSI contours;
% for more information on its structure see the Matlab command 
% 'help contourc'

% intialize the two variables
eddy_lim=[];
largest_curve=[];
% set the maximum number of elements in C
limit = size(C,2);
i = 1; % begin two counters
ii=1;
% rearrange all the contours in C to the structure array 'isolines';
% each element of isolines contains all the vertices of a given contour 
% level of PSI
while(i < limit)
    npoints = C(2,i);
    nexti = i+npoints+1;
    isolines(ii).x = C(1,i+1:i+npoints); % vertex lon's
    isolines(ii).y = C(2,i+1:i+npoints); % vertex lat's
    isolines(ii).l = max(C(2,i+1:i+npoints)); % max lat of a curve
    ii=ii+1;
    i=nexti;
end
% sort the contours accroding to their maximum latitude; this way the first
% closed contour across which velocity increases will also be the largest
% one (it's the one which extend further north).
[unused, order] = sort([isolines(:).l],'descend');
isolines=isolines(order);
% restart the counter and initialize the two flags
i=1;
closed_indx=0; % set to 1 when eddy shape is found
largest_indx=0; % set to 1 when largest closed contour is found 
% inspect all isolines until the eddy shape is determined
% (closed_indx=1 stops the loop)
while closed_indx==0 && i<=length(isolines)
    xdata=isolines(i).x; % vertex lon's
    ydata=isolines(i).y; % vertex lat's
    % conditions to have the largest closed contour around the center
    % (isolines already sorted by maximum latitude)
    % 1) closed contours
    % 2) detected eddy center inside the polygon
    if xdata(1)==xdata(end) && ydata(1)==ydata(end) && ...
            inpolygon(km_ci,km_cj,xdata,ydata)
        % find the contour extremes
        Nj=max(ydata); Ni=max(xdata(ydata==Nj));
        Sj=min(ydata); Si=min(xdata(ydata==Sj));
        Ei=max(xdata); Ej=min(ydata(xdata==Ei));
        Wi=min(xdata); Wj=max(ydata(xdata==Wi));
        % check if velocity across the contour increases
        dir=['N';'S';'E';'W'];
        pts_I=[Ni,Si,Ei,Wi];
        pts_J=[Nj,Sj,Ej,Wj];
        % inspect one extreme at the time (faster)
        ii=1; % counter
        smaller_vel=0; % flag to stop the loop (1 if velocity decreases 
                       % across on of the extremes)
        while ii<=length(dir) && smaller_vel==0
            % interpolate velocity across the extreme
            in_out_vel=interpolate_limits(pts_I(ii),pts_J(ii), ...
                km_i,km_j,dkm,vel,dir(ii));
            % change the flag value if velocity decreases 
            if in_out_vel(1,1)>in_out_vel(2,1)
                smaller_vel=1;
            end
            ii=ii+1; % increase the counter
        end        
        % only if velocity increases across all four extremes the closed 
        % contour is saved as eddy shape
        if smaller_vel==0
            eddy_lim=[xdata;ydata];
            closed_indx=1;
        end
        % largest closed conotur is saved as well
        if largest_indx==0
            largest_curve=[xdata;ydata];
            largest_indx=1;
        end
    end
    i = i+1; % increase the counter
end
% in case velocity doesn't increase across the closed contour, eddy shape
% is defined simply as the largest closd contour
if isempty(eddy_lim) && ~isempty(largest_curve)
    eddy_lim=largest_curve;
    large=1;
end
% eddy_lim is empty if no closed contours of psi exist. In that case eddy
% shape is assumed to be circular with radius 'a-1'. See eddy_dim.m for
% further details.
