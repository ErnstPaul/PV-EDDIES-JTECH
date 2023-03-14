function [Tx,Ty]=ra_windstr(u,v)
% Tau=Rho*Cd*(speed)^2; Tx=Rho*Cd*Speed*u; Ty=Rho*Cd*Speed*v
%===========================================================%
% RA_WINDSTR  $Id: ra_windstr.m, 2014/10/29 $
%          Copyright (C) CORAL-IITKGP, Ramkrushn S. Patel 2014.
%
% AUTHOR: 
% Ramkrushn S. Patel (ramkrushn.scrv89@gmail.com)
% Roll No: 13CL60R05
% Place: IIT Kharagpur.
% This is a part of M. Tech project, under the supervision of DR. ARUN CHAKRABORTY
%===========================================================%
%
% USAGE: [Tx, Ty]=ra_windstr(u,v)
%  
% DESCRIPTION:  Function to compute wind stress from wind field data Based on Gill, 1982 
% Formula and a non-linear Cd based on Large and Pond (1981), modified for low wind 
% speeds (Trenberth et al., 1990)
% 
% INPUTS: 
% u = Zonal wind component [m/s], must be 2D
% v = Meridional wind component [m/s], must be 2D
%
% OUTPUT: 
% Tx = Zonal wind stress [N/m^2]
% Ty = Meridional wind stress [N/m^2]
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
%
% ACKNOWLEDGMENT:
% Author is eternally grateful to MathWorks for providing in built functions. 
% ***********************************************************************************************%
% Checking for input arguments
if (nargin ~= 2)
    error('ra_windstr.m: requires at least two input arguments')
end
% Checking size of U and V
if (size(u)~=size(v))
    error('ra_windstr.m: SIZE of both wind components must be SAME')
end
% Check output
if (nargout ~= 2)
    error('ra_windstr.m: Output arguments are not sufficient'); end
% Defining Constant 
roh=1.2; % kg/m^3, air density
% Computation of Wind Stresses
[lt, ln]=size(u);
Tx=NaN(lt, ln); 
Ty=NaN(lt, ln);
for ii=1:lt
    for jj=1:ln
        U=sqrt(u(ii, jj)^2+v(ii, jj)^2); % Wind speed
        if (U <= 1)
            Cd=0.00218; 
        else if (U > 1 && U <= 3)
                Cd=(0.62+1.56/U)*0.001;
            else if (U > 3 || U < 10)
                    Cd=0.00114; 
                else 
                    Cd=(0.49+0.065*U)*0.001;
                end % endelseif1
            end % endelseif2
        end% endif
        Tx(ii, jj)=Cd*roh*U*u(ii, jj); % kg/m^3*m/s*m/s= N/m^2
        Ty(ii, jj)=Cd*roh*U*v(ii, jj);
    end% endfor : jj
end% endfor: ii