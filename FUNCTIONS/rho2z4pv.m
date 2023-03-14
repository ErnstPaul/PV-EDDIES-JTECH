function [ Zrho ] = rho2z4pv( kmaxrhoref,rhoref,zref,rhoin )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if ~isnan(rhoin)
k0=find((rhoref-rhoin).^2==min((rhoref-rhoin).^2));
if (k0-1)*(k0-kmaxrhoref)<0
    if rhoin>rhoref(k0)
        Zrho=zref(k0)+(rhoin-rhoref(k0))*(zref(k0+1)-zref(k0))/(rhoref(k0+1)-rhoref(k0));
    else
        Zrho=zref(k0)+(rhoin-rhoref(k0))*(zref(k0-1)-zref(k0))/(rhoref(k0-1)-rhoref(k0));
    end
elseif k0==1
    if rhoin>rhoref(k0)
        Zrho=zref(k0)+(rhoin-rhoref(k0))*(zref(k0+1)-zref(k0))/(rhoref(k0+1)-rhoref(k0));
    else % maybe better to artificially exagerate surface gradient for outcropping ? (use 0.05)
        Zrho=zref(k0)+(rhoin-rhoref(k0))*(zref(k0+1)-zref(k0))/(rhoref(k0+1)-rhoref(k0));%/0.05
    end
else
    if rhoin>rhoref(k0)
        Zrho=zref(k0)+(rhoin-rhoref(k0))*(zref(k0-1)-zref(k0))/(rhoref(k0-1)-rhoref(k0));
    else
        Zrho=zref(k0)+(rhoin-rhoref(k0))*(zref(k0-1)-zref(k0))/(rhoref(k0-1)-rhoref(k0));
    end
end
else
    Zrho=nan;
end
%
end
