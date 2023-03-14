function [matching_contour,area_intersection] = cori_matching_contours(contours,Area,indx,indy,dt,next_map_eddies,tot_dist)

%-----------------------------------------
% Cherche s'il y a une superpositon des contours
% Cherche jusqu'à t+DT cartes
%
% Sorties:
% matching_contour: l'indice des contours qui correspondent, ou [] si rien
% area_intersection: l'aire des intersections, ou [] si rien
%
% Entrées:
% contours = CEs ou AEs (cellules Neddy,Nmaps,2)
% Area = matrice d'aires
% indx = indice du tourbillon testé (ligne)
% indy = indice du numéro de la carte (colonne)
% next_map_eddies = indices des eddies candidats à indy+dt
%
% Juin 2016, Cori Pegliasco
%
%------------------------------------------

% Sort distance in ascending mod and take id position in variable pos
[~, pos] = sort(tot_dist);
% Select the nearest 3 (or less) detected eddies
M = min([3 length(tot_dist)]);

% contour position of the studied eddy
Xcon0 = contours{indx,indy,1};
Ycon0 = contours{indx,indy,2};

if length(pos)>=3
    kk =3 ;
else
    kk = length(pos);
end

Xcon_cell = [];
Ycon_cell = [];
%plot(double(Xcon0),double(Ycon0),'r','linewidth',2)
%hold on
% Contour position of these detected eddies
for i=1:kk
    Xcon_cell{i} = contours{next_map_eddies(pos(i)),indy+dt,1};
    Ycon_cell{i} = contours{next_map_eddies(pos(i)),indy+dt,2};
 %   plot(Xcon_cell{i},Ycon_cell{i})
end

%xlim([18 32])
%ylim([29 36])

% Allow memory to number of point intersected point in x direction
t = zeros(M,1);
area = NaN(M,1);
area2 = NaN(M,1);

% Test each detected eddy comparing contour and storing
% number of point in intersected area (in x direction)
for k=1:M %
    
    [a,b] = poly_intersect2(double(Xcon0),double(Ycon0),double(Xcon_cell{k}),double(Ycon_cell{k}));
    if isempty(a)==0
        t(k)=1;
    end
    area2(k)=calcul_aire_rayon(a,b);
    area(k)=area2(k)*100/double(Area(next_map_eddies(pos(k)),indy+1));
end


if sum(t)==0 % no intersections
    
    matching_contour = [];
    area_intersection = [];
    
elseif sum(t)~=0
    % at least one contour intersection
    
    ind = t~=0;                          % Find indices of eddy(ies) with overlaying contour
    matching_contour = next_map_eddies(pos(ind));    % Position of this eddy(ies) in sorted tables
    area_intersection = area(ind);
end


end