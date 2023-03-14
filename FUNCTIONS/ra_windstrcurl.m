function [curlZ]=ra_windstrcurl(lat,lon,u,v)
% curlZ = dTy/dx - dTx/dy; solved by finite difference method
%===========================================================%
% RA_WINDSTRCURL  $Id: ra_windstrcurl.m, 2015/01/15 $
%          Copyright (C) CORAL-IITKGP, Ramkrushn S. Patel 2014.
%  
% AUTHOR: 
% Ramkrushn S. Patel (ramkrushn.scrv89@gmail.com)
% Roll No: 13CL60R05
% Place: IIT Kharagpur.
% This is a part of M. Tech project, under the supervision of DR. ARUN CHAKRABORTY
%===========================================================%
%
% USAGE: curlz=ra_windstrcurl(lat,lon,u,v)
%  
% PREREQUISITE:
% ra_windstr.m written by same author
%
% DESCRIPTION:  Function to compute wind stress curl from wind field data Based on 
% NRSC (2013)
% 
% INPUTS:
% lat = Latitude vector [deg.]
% lon = Longitude Vector [deg.]
% u = Zonal wind component [m/s], must be 2D
% v = Meridional wind component [m/s], must be 2D
%
% OUTPUT: 
% curlZ = Wind stress curl [N/m^3]
% 
% DISCLAIMER: 
% Albeit this function is designed only for academic purpose, it can be implemented in 
% research. Nonetheless, author does not guarantee the accuracy.
% 
% REFERENCE:
% A.E. Gill, 1982, �Atmosphere-Ocean Dynamics�, Academy Press, Vol. 30.
% W. G. Large & S. Pond., 1981,�Open Ocean Measurements in Moderate to Strong Winds�, 
% J. Physical Oceanography, Vol. 11, pp. 324 - 336.
% K.E. Trenberth, W.G. Large & J.G. Olson, 1990, �The Mean Annual Cycle in Global Ocean 
% Wind Stress�, J.Physical Oceanography, Vol. 20, pp. 1742 � 1760.
% NRSC, 2013, "OSCAT Wind stress and Wind stress curl products", Ocean Sciences Group,
% Earth and Climate Science Area, Hyderabad, India.
%
% ACKNOWLEDGMENT:
% Author is grateful to MathWorks for developing in built functions. 
% ***********************************************************************************************%
% To ensure that function receive ample input arguments
if (nargin ~= 4) 
    error(upper('ra_windstrcurl.m : Input arguments are not sufficient or abundance'))
end
% degree to radian
rad=pi/180;
% wind stresses computation
[Tx,Ty]=ra_windstr(u, v);
% computation of curl
[lt, ln]=size(u);
a=diff(lat);
aa=NaN*ones(length(a)-1,1);
for ii=1:length(a)-1
    if (a(ii) == a(ii+1))
        aa(ii)=a(ii);
    else
        error('Latitude difference is not consistance')
    end % endif
    dlat=mean(aa);
end % endfor
clear ii
deltay=dlat*111176;
curlZ=NaN(lt, ln);
long=NaN(lt, ln);
for ii=1:lt
   for jj=1:ln
           long(ii, jj)=lon(jj)*111176*cos(lat(ii)*rad);
        %long(i,j)=lon(j)*6378137*rad*cos(lat(i)*rad); 
        % [m] earth radious in meters= 6,378,137.0 m.. from wikipedia.
   end % endfor
end % endfor
clear ii jj
% Centeral difference method in x and y 
for ii=2:lt-1
    for jj=2:ln-1
        curlZ(ii, jj)=(Ty(ii, jj+1)-Ty(ii, jj-1))/(2*(long(ii, jj+1)-long(ii, jj-1))) - ...  
            (Tx(ii+1, jj)-Tx(ii-1, jj))/(2*deltay) ;
    end % endfor
end % endfor
clear ii jj
% Forward difference method in x and y 
for jj=1:ln-1
    curlZ(1, jj)=(Ty(1, jj+1)-Ty(1, jj))/(long(1, jj+1)-long(1, jj)) - ...
        (Tx(2, jj)-Tx(1, jj))/deltay ; 
end 
for ii=1:lt-1
    curlZ(ii, 1)=(Ty(ii, 2)-Ty(ii, 1))/(long(ii, 2)-long(ii, 1)) - ...
        (Tx(ii, 2)-Tx(ii, 1))/deltay ;
end
clear ii jj
curlZ(1, ln)=curlZ(1, ln-1);
% Backward difference method in x and y
for ii=2:lt
    curlZ(ii, ln)=(Ty(ii, ln)-Ty(ii, ln-1))/(long(ii, ln)-long(ii, ln-1)) - ...
        (Tx(ii, ln)-Tx(ii-1, ln))/deltay ; 
end
for jj=2:ln-1
    curlZ(lt, jj)=(Ty(lt, jj)-Ty(lt, jj-1))/(long(lt, jj)-long(lt, jj-1)) - ...
        (Tx(lt, jj)-Tx(lt-1, jj))/deltay ;
end
clear ii jj
curlZ(lt, 1)=curlZ(lt, ln-1);