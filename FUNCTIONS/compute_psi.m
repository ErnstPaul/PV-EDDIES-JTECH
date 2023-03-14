function psi=compute_psi(u,v,km_di,km_dj)
% compute_psi.m
%  compute_psi(u,v,km_dlon,km_dlat) computes the streamfunction (PSI) field
%  by spatially integrating the u- and v-component of velocity within the 
%  area around the detected eddy center defined in eddy_dim.m.
%  The underlaying assumption is that in the presence of an eddy the
%  the velocity field is characterized by weak divergence, so that contours
%  of PSI are tangential to the velocity vectors.
%  To reduce the error associated with this assumption PSI is computed
%  integrating u and v along two different paths, and the two fields are
%  then averaged.
%
%  - u and v are NxM matrices of the two component of velocity;
%  - km_di and km_dj are NxM matrices of the longitudinal and latitudinal 
%      spacing in km between grid points. (They can vary along i and j)
%
%  OUTPUT:
%  - psi is the NxM matrix of the streamfunction field 
%
%  See also Appendix in Nencioli et al. (2009) for further details.
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

% Indices for domain size
lx=size(u,1);
ly=size(u,2);
% adjust km_di and km_dj for trapz integration
% (add first row/column of zeros)
km_di=[zeros(lx,1) km_di];
km_dj=[zeros(1,ly); km_dj];
% itegrate first row of v along longitude (first term of eq.A2)
cx=cumtrapz(v(1,:)).*km_di(1,:); % trapezoidal sum
% integrate first column of u along latitude (second term of eq.A3)
cy=cumtrapz(u(:,1)).*km_dj(:,1);
% expand the vectors into matrices to compute PSI
cx=cx(ones(lx,1),:);
cy=cy(:,ones(1,ly));
% compute streamfunction -----------------------
% PSI from integrating v firts and then u
psi_xy=(-cx + cumtrapz(u).*km_dj); %(eq. A2)
% PSI from integrating u first and then v 
psi_yx=(-cumtrapz(v,2).*km_di + cy); %(eq. A3)
% final PSI as average between the two
psi=(psi_xy+psi_yx)/2;
%---------------------------------------------

%psi=((cumtrapz(v,2).*km_dlon-cy)-(cumtrapz(u)*km_dlat-cx))/2;
