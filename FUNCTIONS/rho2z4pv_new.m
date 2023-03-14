function [ Zrho ] = rho2z4pv( kmaxrhoref,rhoref,zref,rhoin )
%
%   Calculate level associated with density rhoin 
%   in the reference profile
%   rhoref
%
epsrho=1e-11;% minimum drho for interpolation division
%
if ~isnan(rhoin)

    if rhoin<rhoref(1)
        k0=1;
        % extrapolate value : maybe better to artificially exagerate surface gradient for outcropping ? (use 1o/oo over 50 m  ...)
        Zrho=zref(k0)+(rhoin-rhoref(k0))*50.;%
        %Zrho=zref(k0)+(rhoin-rhoref(k0))*(zref(k0-1)-zref(k0))/(rhoref(k0-1)-rhoref(k0));%/0.05
    elseif rhoin>rhoref(kmaxrhoref)
        k0=kmaxrhoref;
        % extrapolate value : maybe better to artificially exagerate gradient for outcropping ? (use 1o/oo over 50 m  ...)
        Zrho=zref(k0)+(rhoin-rhoref(k0))*50.;%
        %Zrho=zref(k0)+(rhoin-rhoref(k0))*(zref(k0-1)-zref(k0))/(rhoref(k0-1)-rhoref(k0));%/0.05
    else
        k0=min(find((rhoref-rhoin).^2==min((rhoref-rhoin).^2)));
        if rhoin==rhoref(k0)
          Zrho=zref(k0);
        elseif rhoin>rhoref(k0)
          Zrho=zref(k0)+(rhoin-rhoref(k0))*(zref(k0+1)-zref(k0))/max(rhoref(k0+1)-rhoref(k0),epsrho);
        else
          Zrho=zref(k0)+(rhoin-rhoref(k0))*(zref(k0)-zref(k0-1))/max(rhoref(k0)-rhoref(k0-1),epsrho);
        end
    end
%
else
    Zrho=nan;
end
%
end
