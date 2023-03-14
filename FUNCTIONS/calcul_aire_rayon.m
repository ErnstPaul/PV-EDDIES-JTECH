function [Aire,Rayon] = calcul_aire_rayon(Xcon,Ycon);

% calcule l'aire et le rayon equivalent d'une surface sur la terre
%
% inputs:  Xcon,Ycon = coordonnees en longitude/latitude du contour de la surface (obtenues par convhull
%                         ou autre...)
% outputs:   Aire en km^2, rayon en km


dy = 111.12;  % 1º de latitude = 111.12 km 

X = Xcon.*cos(Ycon*pi/180)*dy;
Y = Ycon.*dy;

Aire = polyarea(X,Y);
Rayon = sqrt(Aire/pi);
