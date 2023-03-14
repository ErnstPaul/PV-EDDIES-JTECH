function [Anticyclonic_Cell,Cyclonic_Cell,Ncyclo,Nanti] = Threshold_eddy_edges(FIELD,STATS,X,Y,FACTOR,EXTREMES,minRad,minAmp,minRatio)
%Returns the cells as in EddyExtractionMaster but for thresholds
%=========================
%----------------------------------------
% initialise cells
Cyclonic_Cell=cell(length(EXTREMES),12);
Anticyclonic_Cell=cell(length(EXTREMES),12);
Ncyclo = 0;
Nanti = 0;

%---------------------------------------------
% Define contour, including potential centers
LOW = abs(FIELD);
LOW(isnan(FIELD)) = 0;

% Create regular grid
Xgrid = zeros(length(X),length(Y));
Ygrid = zeros(length(X),length(Y));
for p = 1:length(Y)
    Xgrid(:,p) = X(:,1);
end
for p = 1:length(X)
    Ygrid(p,:) = Y(:,1);
end

%The contour factor is calculated here
K = STATS(1,2)*FACTOR;

%get contours
HF = figure('visible','off');
CS = contour(Xgrid,Ygrid,LOW,[K K]);
close(HF)
X = Xgrid;
Y = Ygrid;


% Initialization
j = 1; % first coordinates of the contour scan
k = 1; % first contour
n_min = 10; %minimum 10 grid points (4 inside)

% scan each FIELD contour
while j < size(CS,2)

    n = CS(2,j); % number of coordinates for the contour(j)
    xd = CS(1,j+1:j+n); % X values serie for the contour(j) coordinates
    yd = CS(2,j+1:j+n); % Y values serie for the contour(j) coordinates

    % validate only bigger contour
    if n >= n_min

        % is a closed contour
        if xd(end) == xd(1) && yd(end) ==yd(1)
            % make a mask of the contour
            in = InPolygon(X,Y,xd,yd);
            Lm = FIELD;
            Lm(~in) = NaN;
            %find extremes in this contour
            toCheck = find(ismember(EXTREMES,find(double(in))));
            %if we got nothing here
            if (isempty(toCheck))
                j = j + n + 1; % series of coordinates of the next contour
                continue;
            end
            %select only contours with one extreme here
            if (length(abs(Lm(EXTREMES(toCheck)))) > 1)
                j = j + n + 1; % series of coordinates of the next contour
                continue;
            end
            %In case this is screwed up
            exVal = nanmax(abs(Lm(EXTREMES(toCheck))), [], 'all');
            if (isnan(exVal))
                j = j + n + 1; % series of coordinates of the next contour
                continue;
            end
            exInd = find(abs(Lm)==exVal);
            exInd = exInd(1);

            %Check for edgefield as well to ensure eddy shape
            LC = Lm(abs(Lm)==max(abs(Lm(:))));
            xLmax = X(find(Lm==LC(1),1));
            yLmax = Y(find(Lm==LC(1),1));

            %get X and Y indices
            [xc,yc] = ind2sub(size(Lm),exInd);
            exVal = FIELD(xc,yc);

            %check for amplitude
            if (abs(exVal) >= minAmp)

                % calculate area and radius
                [b,a] = calcul_aire_rayon(xd,yd);

                %check for radius
                if (a >= minRad)

                    %compute the ratio of axes
                    geom = polygeom(xd,yd);
                    centroidX = geom(2); centerX = Xgrid(xc);
                    centroidY = geom(3); centerY = Ygrid(yc);
                    distMat = zeros(3,length(xd));
                    for lenCont = 1:length(xd)
                        distMat(1,lenCont) = sqrt((xd(lenCont)-centroidX)^2 + (yd(lenCont)-centroidY)^2);
                        distMat(2,lenCont) = sqrt((xd(lenCont)-centerX)^2 + (yd(lenCont)-centerY)^2);
                        distMat(3,lenCont) = sqrt((xd(lenCont)-xLmax)^2 + (yd(lenCont)-yLmax)^2);
                    end
                    minorAxisOID = min(distMat(1,:));
                    majorAxisOID = max(distMat(1,:));
                    ratioAxesOID = minorAxisOID/majorAxisOID;
                    minorAxis = min(distMat(2,:));
                    majorAxis = max(distMat(2,:));
                    ratioAxes = minorAxis/majorAxis;
                    minorAxisED = min(distMat(3,:));
                    majorAxisED = max(distMat(3,:));
                    ratioAxesED = minorAxisED/majorAxisED;

                    %check for ratios between centers and centroid
                    if ((ratioAxesOID>minRatio) && (ratioAxes>minRatio/2) && (ratioAxesED>minRatio))

                        %cyclone vs anticyclone determination: center greater than the mean+K or less
%                         disp(['cent ' num2str(exVal) ' mean ' num2str((STATS(1,1)+K))])
                        if (exVal > (STATS(1,1)+K))

                            Nanti = Nanti+1;

                            % Fill values of X and Y center
                            Anticyclonic_Cell{Nanti,1}=centerX;
                            Anticyclonic_Cell{Nanti,2}=centerY;

                            % Fill values of X and Y  geoïd
                            Anticyclonic_Cell{Nanti,3}=single(geom(2));
                            Anticyclonic_Cell{Nanti,4}=single(geom(3));

                            % Fill tables of X and Y contours
                            Anticyclonic_Cell{Nanti,5}=single(xd);
                            Anticyclonic_Cell{Nanti,6}=single(yd);

                            % Geometric Area and Radius
                            Anticyclonic_Cell{Nanti,7}=single(a);
                            Anticyclonic_Cell{Nanti,8}=single(b);

                            % Fill values of Amplitude as fin = h(center)-h(contour)
                            Anticyclonic_Cell{Nanti,9}=single(exVal-K);

                            % Fill ratio values
                            Anticyclonic_Cell{Nanti,11}=single(ratioAxesOID);
                            Anticyclonic_Cell{Nanti,12}=single(ratioAxes);

                        else %is cyclone
                            Ncyclo = Ncyclo+1;

                            % Fill values of X and Y center
                            Cyclonic_Cell{Ncyclo,1}=centerX;
                            Cyclonic_Cell{Ncyclo,2}=centerY;

                            % Fill values of X and Y  geoïd
                            Cyclonic_Cell{Ncyclo,3}=single(geom(2));
                            Cyclonic_Cell{Ncyclo,4}=single(geom(3));

                            % Fill tables of X and Y contours
                            Cyclonic_Cell{Ncyclo,5}=single(xd);
                            Cyclonic_Cell{Ncyclo,6}=single(yd);

                            % Geometric Area and Radius
                            Cyclonic_Cell{Ncyclo,7}=single(a);
                            Cyclonic_Cell{Ncyclo,8}=single(b);

                            % Fill values of Amplitude as fin = h(center)-h(contour)
                            Cyclonic_Cell{Ncyclo,9}=single(exVal-K);

                            %Fill ratios
                            Cyclonic_Cell{Ncyclo,11}=single(ratioAxesOID);
                            Cyclonic_Cell{Ncyclo,12}=single(ratioAxes);
                        end %endAEorCEIF
                    end %endratioIF
                end %endminradIF
            end %endminampIF
        end %endclosedcontourIF
    end %endbiggercontourIF
    % increment the counter
    j = j + n + 1; % series of coordinates of the next contour
end

%Truncate 
Anticyclonic_Cell(Nanti+1:end,:) = [];
Cyclonic_Cell(Ncyclo+1:end,:) = [];
