function[xo,yo] = poly_intersect(x1,y1,x2,y2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%    Find xo and yo, the successive points constituting the intersection
%%%    between the 2 polygons (x1,y1) and (x2,y2)
%%%
%%%    Use AREAGK, CORI_ISINTPL, INPOLYGON, CALCUL_AIRE_RAYON, AC_DISTANCE
%%%    June 2016  Cori Pegliasco
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



flag = 1;

% Handle input
if nargin < 4
    error(' Not enough input arguments')
end

x1 = x1(:); y1 = y1(:);
x2 = x2(:); y2 = y2(:);
l1 = length(x1);
l2 = length(x2);

if length(y1)~=l1 | length(y2)~=l2
    error('(X1,Y1) and (X2,Y2) must pair-wise have the same length.')
end

% Check areas and reverse if negative
nn1 = areagk(x1,y1);
if nn1<0, x1 = flipud(x1); y1 = flipud(y1); end
nn2 = areagk(x2,y2);
if (nn2<0 & flag<3) | (nn2>0 & flag==3)
    x2 = flipud(x2); y2 = flipud(y2);
end

% If both polygons are identical ........
if l1==l2
    if all(x1==x2) & all(y1==y2)
        if flag<3, xo = x1; yo = y1; ind = 1:l1;
        else, xo = []; yo = []; ind = []; end
        return
    end
end

% Only 4 points of a contour inside the other -> quick exit

one_in1 = inpolygon(x1,y1,x2,y2);
one_in2 = inpolygon(x2,y2,x1,y1);
ii1 = find(one_in1 == 1);
ii2 = find(one_in2 == 1);

if length([ii1; ii2])<4
    xo = []; yo = [];
    return
end

% Calculate matrix of intersections .....

[is,is2,C] = cori_isintpl(x1,y1,x2,y2);



% cas des tourbillons imbriqués et cas de pas d'intersections
if ~is2 && is
    % 2 cas: petit tourbillon dans grand
    % tourbillons séparés
    area1 = calcul_aire_rayon(x1,y1);
    area2 = calcul_aire_rayon(x2,y2);
    [~,ind_area_min] = min([area1 area2]);
    [~,ind_area_max] = max([area1 area2]);
    
    xx{1} = x1;
    xx{2} = x2;
    yy{1} = y1;
    yy{2} = y2;
    is_inside = [];
    for i=1:length(xx{ind_area_min})
        is_inside(i) = inpolygon(xx{ind_area_min}(i),yy{ind_area_min}(i),xx{ind_area_max},yy{ind_area_max});
    end
    ind_1 = [];
    ind_1 = find(is_inside ==1);
    if isempty(ind_1)==0
        xo=xx{ind_area_min}; yo = yy{ind_area_min};
        return
    end
    xo = [];
    yo = [];
    
    return
elseif ~is2 && ~is   % Quick exit if no intersections ........
    xo = [];
    yo = [];
    
end

% Mark intersections with unique numbers
i1 = find(C);
ni = length(i1);
C(i1) = 1:ni;


% Calculate intersections themselves
[i1,i2,id] = find(C~=0);

if isempty(i1)==1
    xo = [];
    yo = [];
    return
else
    
    xtt = [x1(ii1); x2(ii2)];
    ytt = [y1(ii1); y2(ii2)];
    
    
    
    xbar = sum(xtt)/length(xtt);
    ybar = sum(ytt)/length(ytt);
    angle = atan2 ( (ytt - ybar), (xtt - xbar) );
    [~,pos] = sort(angle);
    xt = xtt(pos);
    yt = ytt(pos);
    
    % supprimer les points doubles
    indx0 = find(diff(xt)==0);
    indy0 = find(diff(yt)==0);
    ind0 = ismember(indx0,indy0);
    xt(indx0(ind0)) = [];
    yt(indx0(ind0)) = [];
    
    % jointer le contour
    xo = [xt; xt(1)];
    yo = [yt; yt(1)];
    
    % figure
    % hold on
    % plot(x1,y1,'rx-')
    % plot(x2,y2,'bx-')
    % for j=1:length(xt)
    %     plot(xt(j),yt(j),'sk')
    %     pause
    % end
    
    
    
end


%     % calculer pour chaque segment les distances entre les intersections pour
%     % garder les distances (pas les différences d'indices) les plus courtes
%     itot1 = [1; sort(i1); length(x1)];
%     itot2 = [1; sort(i2); length(x2)];
%
%     i1 = sort(i1);
%     i2 = sort(i2);
%     dist1 = NaN([1 length(itot1)]);
%     dist2 = NaN([1 length(itot2)]);
%     dist1sum = NaN([1 length(itot1)-1]);
%     dist2sum = NaN([1 length(itot2)-1]);
%
%     for i=1:length(x1)-1 % Calcul de chaque distance
%         dist1(i) = ac_distance(y1(i),x1(i),y1(i+1),x1(i+1));
%     end
%     for i=1:length(x2)-1
%         dist2(i) = ac_distance(y2(i),x2(i),y2(i+1),x2(i+1));
%     end
%
%     for i=1:length(itot1)-1 % Somme des distances entre les intersections
%         dist1sum(i) = sum(dist1(itot1(i):itot1(i+1)-1));
%     end
%     for i=1:length(itot2)-1
%         dist2sum(i) = sum(dist2(itot2(i):itot2(i+1)-1));
%     end
%     % réorganisation des sommes
%     disttot1 = [dist1sum(2:length(itot1)-2) dist1sum(1)+dist1sum(end)];
%     disttot2 = [dist2sum(2:length(itot2)-2) dist2sum(1)+dist2sum(end)];
%
%     itot(1,:) = [i1; length(x1); 1; i1(1)];
%     itot(2,:) = [i2; length(x2); 1; i2(1)];
%     xx{1} = x1;
%     xx{2} = x2;
%     yy{1} = y1;
%     yy{2} = y2;
%
%     [~,posmindist] = min([disttot1; disttot2]); % cherche si la distance min est sur le polygone 1 ou 2
%     xcon = [];
%     ycon = [];
%
%     for i=1:length(itot)-1
%
%         if i<length(itot)-2 % write points between intersections
%             xcon =  [xcon; xx{posmindist(i)}(itot(posmindist(i),i) : itot(posmindist(i),i+1))];
%             ycon =  [ycon; yy{posmindist(i)}(itot(posmindist(i),i) : itot(posmindist(i),i+1))];
%         else % close the contour (from the last intersection to the first one)
%             xcon =  [xcon; xx{posmindist(end)}(itot(posmindist(end),i): itot(posmindist(end),i+1))];
%             ycon =  [ycon; yy{posmindist(end)}(itot(posmindist(end),i): itot(posmindist(end),i+1))];
%         end
%     end
%
%     % supprimer les points doubles
%     % remove similar points
%     indx0 = find(diff(xcon)==0);
%     indy0 = find(diff(ycon)==0);
%     ind0 = ismember(indx0,indy0);
%     xcon(indx0(ind0)) = [];
%     ycon(indx0(ind0)) = [];
%
%     xo = xcon;
%     yo = ycon;
%     return
%     % xo et yo sont les contours de l'intersection (quand il y a crossing)
%








