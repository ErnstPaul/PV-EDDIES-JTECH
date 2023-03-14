function [next_map_eddies,tot_dist,found_centroid_at] = cori_find_eddy_next_map(Xcentroid,Ycentroid,indx,indy,dt,dist_lim,nc)
t = indy;


%-------------------------------------------
%
%  Cherche les tourbillons sur la carte map+dt
%
%  Entrées: 
%  Xcentroid, Ycentroid: matrices ayant Nx lignes et Ny colonnes
%  indx = indice tourbillon (lignes)
%  indy = indice du numéro de la carte (colonnes)
%  dt = pas de temps entre 2 cartes
%  Au temps t on connait la position du centroide du tourbillon étudié
%  Au temps t+dt on cherche la position des centroides des
%  tourbillons à moins de la distance dist_lim du centroid du temps t
%
%  Sortie: Les indices des tourbillons correspondants, et la distance entre
%  les centroides à t+dt et le centroid étudié
%
%  Juin 2016, Cori Pegliasco
%
%-------------------------------------------

% condition: t+dt <= Nbre max de cartes (nc)

if t+dt > nc
    next_map_eddies = [];
    tot_dist = [];
    found_centroid_at = NaN;
    
    return
end


next_map_eddies = [];
tot_dist = [];
next_map_eddies = find(Xcentroid(:,indy+dt)>=Xcentroid(indx,indy)-dist_lim & Xcentroid(:,indy+dt)<=Xcentroid(indx,indy)+dist_lim ...
    & Ycentroid(:,indy+dt)>=Ycentroid(indx,indy)-dist_lim & Ycentroid(:,indy+dt)<=Ycentroid(indx,indy)+dist_lim );

tot_dist = ac_distance(Ycentroid(indx,indy),Xcentroid(indx,indy),Ycentroid(next_map_eddies,indy+dt),Xcentroid(next_map_eddies,indy+dt));

found_centroid_at = dt;

end






