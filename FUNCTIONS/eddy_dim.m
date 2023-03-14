function [lonlat, warn, land, box, large]=eddy_dim(u,v, ...
    C_I,lon,lat,Lonmin,Lonmax,Latmin,Latmax, ...
    fac,rad,a)
% eddy_dim.m
%  eddy_dim(U,V,day,C_I,C_J,lon,lat,Lonmin,Lonmax,Latmin,Latmax,mask,
%  fac,rad,a) computes the shape of the eddy defined by the center C_I,C_J.
%
%  - U and V are the time series of the 2D u and v velocity field;
%  - day is the day of the 2D velocity field to analize;
%  - C_I and C_J are the indices of the eddy center
%  - lon and lat are longitudes and latitudes of points in the domain;
%  - Lonmin, Lonmax, Latmin and Latmax are the extremes of the domain;
%  - mask is the matrix that defines sea (1) and land points (0);
%  - fac is the factor to enlarge the area where the shape is computed;
%  - rad is used to define the area where the shape is computed;
%  - a is one of the parameters from the eddy detection; is used here to 
%    define eddy radius in case no closed contour of PSI is found;
%
%  OUTPUT:
%  - lonlat is the array containing longitudes (first row) and
%    latitudes(second row) of the vertices defining the eddy shape;
%  - warn, land, box and large are flags that give some informations on the
%    procedure used to compute the eddy shape. See 'mod_eddy_shapes.m' for
%    further details.
%
% 'compute_psi' is used to compute the streamfunction field integrating u-
% and v-component f velocity. Check the documentation in 'compute_psi.m' 
% for further details.
%
% 'max_curve' is used to compute eddy shape (defined as the largest closed 
% contour of PSI around the eddy center across which velocity magnitude 
% increases). Check the documentation in 'max_curve.m' for further details.
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

% increase the dimensions of the area where eddy shape is computed
rad=rad*fac;
% initialized the flags
land=0;
warn=0;
box=0;
%-----------------------------------------------------------
% Added option for periodic East-West boundaries 
% (first and last column of the domain are identical)
%grab grid
Xgrid = zeros(length(lon),length(lat));
Ygrid = zeros(length(lon),length(lat));
for p = 1:length(lat)
    Xgrid(:,p) = lon(:,1);
end
for p = 1:length(lon)
    Ygrid(p,:) = lat(:,1);
end
lat = Ygrid; lon = Xgrid;

% center coordinates (C_J and C_I are center indices in the whole domain) 
c_lat=lat(C_I);
c_lon=lon(C_I);
[C_J, C_I] = ind2sub(size(Xgrid),C_I);

% resize coordinate and velocity matrix 
% (making sure not to go outside the domain)
lat=lat(max(C_J-rad,1):min(C_J+rad,size(lat,1)), ...
    max(C_I-rad,1):min(C_I+rad,size(lat,2)));
lon=lon(max(C_J-rad,1):min(C_J+rad,size(lon,1)), ...
    max(C_I-rad,1):min(C_I+rad,size(lon,2)));
v=v(max(C_J-rad,1):min(C_J+rad,size(v,1)), ...
    max(C_I-rad,1):min(C_I+rad,size(v,2)));
u=u(max(C_J-rad,1):min(C_J+rad,size(u,1)), ...
    max(C_I-rad,1):min(C_I+rad,size(u,2)));

% indices of the eddy center in the smaller area
[c_j, c_i]=find(lat==c_lat & lon==c_lon);

% inspect if there are land-points in the area:
% if land-points are present, the area is resized to have only ocean-points
% (rad becomes the shortes distance from the center to land)
[yl, xl]=find(isnan(u));
if ~isempty(xl)
    land=1; % update flag
    % compute new 'rad' --------------------------------------
    % distance in grid points from the center to the land points
    dpts=sqrt((yl-c_j).^2+(xl-c_i).^2);
    % problem if the land point is a vertex of the region
    if floor(min(dpts)/sqrt(2))==min(dpts)/sqrt(2)
        rad2=floor(min(dpts)/sqrt(2))-1;
    else
        rad2=floor(min(dpts)/sqrt(2));
    end
    % rad2 cannot be less than a-1 (bc of fourth constraint)
    if rad2<a-1
        rad2=a-1;
    end
    % update 'rad' (rad2 has to be bigger than rad at previous fac)
    if rad2>rad/fac*(fac-1)
        rad=rad2;
    else
        rad=rad/fac*(fac-1);
    end
    %---------------------------------------------------------
    % resize coordinates and velocities
    lat=lat(max(c_j-rad,1):min(c_j+rad,size(lat,1)), ...
        max(c_i-rad,1):min(c_i+rad,size(lat,2)));
    lon=lon(max(c_j-rad,1):min(c_j+rad,size(lon,1)), ...
        max(c_i-rad,1):min(c_i+rad,size(lon,2)));    
    v=v(max(c_j-rad,1):min(c_j+rad,size(v,1)), ...
        max(c_i-rad,1):min(c_i+rad,size(v,2)));
    u=u(max(c_j-rad,1):min(c_j+rad,size(u,1)), ...
        max(c_i-rad,1):min(c_i+rad,size(u,2)));
    % indices of the eddy center in the new smaller area
    [c_j c_i]=find(lat==c_lat & lon==c_lon);
end
% compute velocity magnitude
vel=sqrt(u.^2+v.^2);

% convert lon and lat into km distance matrices 
% (needed to compute the streamfunction field)

for i=1:size(lon,1)
    km_di(i,:)=m_lldist(lon(i,:),lat(i,:));
    km_i(i,:)=[0 cumsum(km_di(i,:))];
end
for i=1:size(lat,2)
    km_dj(:,i)=m_lldist(lon(:,i),lat(:,i));
    km_j(:,i)=[0; cumsum(km_dj(:,i))];
end

% compute grid km spacing -----------
% (used for psi computation and 
%  for conversion of eddy boundary from km to lat lon)
%
% conversion factor from km to deg lat (equal lat spacing)
km_dlat=m_lldist([lon(1,1) lon(1,1)],[lat(1,1) lat(1,2)]);
km2lat=(lat(1,2)-lat(1,1))/km_dlat;
% conversion factor from km to deg lon (varying lon spacing with lat)
for i=1:size(lon,1)
    km_dlon(i,1)=m_lldist([lon(i,1) lon(i,2)],[lat(i,1) lat(i,2)]);
end
dlon=diff(lon,1,2);
dlon=dlon(:,1);
km2lon=dlon./km_dlon;
%------------------------------------

% compute streamfunction field
psi=compute_psi(u,v,km_di,km_dj);
% eddy center position in the small area (kilometric coordinates)
km_cj=km_j(c_j,c_i);
km_ci=km_i(c_j,c_i);
% compute eddy shape
[eddy_lim, large]=max_curve(km_i,km_j,psi,km_ci,km_cj,vel);

% in case there is no streamfunction curve closed around the center the
% eddy dimension is set to a-1 points radius
if isempty(eddy_lim)
    warn=1; % update flag
    R=min([km_i(1,1+(a-1))-km_i(1,1) km_i(end,1+(a-1))-km_i(end,1) ...
        km_j(1+(a-1),1)-km_j(1,1) km_j(end,1)-km_j(end-(a-1),1)]);
    theta = 0:pi/20:2*pi;
    lon = (R*cos(theta)+km_ci);
    lat = (R*sin(theta)+km_cj);
    eddy_lim=[lon;lat];
    % in case eddy shape is too close to the area boundary, R is further
    % reduced to 1
    if ~isempty(find(eddy_lim<0 , 1))
        R=min([km_i(1,2)-km_i(1,1) km_i(end,2)-km_i(end,1) ...
        km_j(2,1)-km_j(1,1) km_j(end,1)-km_j(end-1,1)]);
        theta = 0:pi/20:2*pi;
        lon = (R*cos(theta)+km_ci);
        lat = (R*sin(theta)+km_cj);
        eddy_lim=[lon;lat];
    end
end

% parameter that defines if the computed eddy shape is too close to the 
% area boundaries.
Xgrid = zeros(length(lon),length(lat));
Ygrid = zeros(length(lon),length(lat));
for p = 1:length(lat)
    Xgrid(:,p) = lon(1,:);
end
for p = 1:length(lon)
    Ygrid(p,:) = lat(1,:);
end
lat = Ygrid; lon = Xgrid;
lim_bound=1.5 * ...      % distance is set in km
    m_lldist([lon(c_j,1) lon(c_j,2)],[lat(c_j,1) lat(c_j,1)]); 
% if the eddy shape is less than 'lim_bound' from teh area boundaries, and
% the area boundaries does not coincides with the domain boundaries, then
% the flag box is set to 1 so that eddy shape will be computed in a larger
% area (see also mod_eddy_shapes.m)
if (min(eddy_lim(1,:))-min(min(km_i)) < lim_bound || ...
        max(max(km_i))-max(eddy_lim(1,:)) < lim_bound || ...
        min(eddy_lim(2,:))-min(min(km_j)) < lim_bound || ...
        max(max(km_j))-max(eddy_lim(2,:)) < lim_bound) && ...
        (min(min(lon))~=Lonmin && ...
        max(max(lon))~=Lonmax && ...
        min(min(lat))~=Latmin && ...
        max(max(lat))~=Latmax)
    box=1;   
end

% convert eddy_lim from km to geographic coordinates -----------
% initialize lonlat
lonlat=zeros(size(eddy_lim));
% convert point by point
for j=1:length(eddy_lim(1,:))
	% Find i and j index of closest grid point to eddy boundary point
	aa=find(km_i<=eddy_lim(1,j) & km_j<=eddy_lim(2,j));
	minD=sqrt((km_i(aa)-eddy_lim(1,j)).^2+(km_j(aa)-eddy_lim(2,j)).^2);
	[J, I]=ind2sub(size(km_i),aa(minD==min(minD)));
	% Closest point position in fraction of index
	if I==size(km_i,2) % To fix case where contour is out of domain 
			   % (Matlab bug?)
		delta_i=1;
		I=I-1;
	else
		delta_i=(eddy_lim(1,j)-km_i(J,I))/(km_i(J,I+1)-km_i(J,I));
	end
	if J==size(km_j,1) % To fix case where contour is out of domain 
		           % (Matlab bug?)
		delta_j=1;
		J=J-1;
	else
		delta_j=(eddy_lim(2,j)-km_j(J,I))/(km_j(J+1,I)-km_j(J,I));
	end
	% Lon and Lat coordinates of eddy boundary point
	% (bi-linear interpolator)
	lonlat(1,j)=lon(J,I)*(1-delta_i)*(1-delta_j)+ ...
		lon(J,I+1)*delta_i*(1-delta_j)+ ...
		lon(J+1,I)*(1-delta_i)*delta_j+ ...
		lon(J+1,I+1)*delta_i*delta_j;
	lonlat(2,j)=lat(J,I)*(1-delta_i)*(1-delta_j)+ ...
		lat(J,I+1)*delta_i*(1-delta_j)+ ...
		lat(J+1,I)*(1-delta_i)*delta_j+ ...
		lat(J+1,I+1)*delta_i*delta_j;
	% i and j coordinates of eddy boundary point
	lonlat(3,j)=I+delta_i-c_i+C_I;
	lonlat(4,j)=J+delta_j-c_j+C_J;
end
%---------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% uncomment to plot the shapes of the eddies detected in the domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if large==1 || warn==1
%     clr='r';
% else
%     clr='k';
% end
% figure
% quiver(km_i,km_j,u,v,'k')
% hold on
% plot(km_ci,km_cj,'k*')
% contour(km_i,km_j,psi,100)
% line(eddy_lim(1,:),eddy_lim(2,:),'linewi',2,'color',clr);
% hold off
% 
% if ~isempty(lonlat)
%     figure
%     quiver(lon,lat,u,v,'k')
%     hold on
%     plot(c_lon,c_lat,'r*')
%     contour(lon,lat,psi,100)
%     line(lonlat(1,:),lonlat(2,:),'linewi',2,'color',clr);
%     hold off
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
