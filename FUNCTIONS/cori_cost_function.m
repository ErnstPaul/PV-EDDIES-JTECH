function [min_CF] = cori_cost_function(matching_contour,indy,indx,Radius,Amplitude,EKE,Xcenter,Ycenter,...
    mean_delt_Rad,mean_delt_Ampl,mean_delt_EKE,mean_delt_dist,std_delt_Rad,std_delt_Ampl,std_delt_EKE,std_delt_dist,DT)

%------------------------------------
%
% Minimise la fonction de coût définie dans Pegliasco et al [2015]
% Rajout de la distance entre la position des centres consécutifs
%
% Entrées:
% matching_contour: contours en intersection avec l'eddy considéré
% indx,indy: indices des tourbillons (lignes) et du n° de carte (colonnes)
% Radius, Amplitude, EKE, Xcenter, Ycenter: les matrices organisées par Amplitude décroissante
% mean_xx et std_xx: statistiques sur les trajectoires
%
% Sortie: indice de la position du minimum de la fonction de cout
%
% Juin 2016, Cori Pegliasco
%
%------------------------------------

cf_radius=(((Radius(matching_contour,indy+DT)-Radius(indx,indy))- mean_delt_Rad)./std_delt_Rad).^2;
cf_eke=(((EKE(matching_contour,indy+DT)-EKE(indx,indy))- mean_delt_EKE)./std_delt_EKE).^2;
cf_amp=(((Amplitude(matching_contour,indy+DT)-Amplitude(indx,indy))- mean_delt_Ampl)./std_delt_Ampl).^2;
cf_dist = ((ac_distance(Ycenter(indx,indy),Xcenter(indx,indy),Ycenter(matching_contour,indy+DT),Xcenter(matching_contour,indy+DT))- mean_delt_dist)./std_delt_dist).^2;
cf =  cf_radius + cf_eke + cf_amp + cf_dist;
[~,min_CF] = min(cf);



end
