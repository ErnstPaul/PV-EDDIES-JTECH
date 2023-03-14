%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calcul du parametre de coriolis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f = coriolis(lat)

f = 2*2*pi*sin(lat*pi/180)/(24*3600);