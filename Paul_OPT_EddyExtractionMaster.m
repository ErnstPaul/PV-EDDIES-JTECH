%% Paul_OPT_EddyExtractionMaster.m (version 2.2)
%Author: Paul Ernst
%Date Created: 4/27/2022
%Date of Last Update: 1/12/2023
%Update History:
%PE 1/12/2023 - Cleared in final revision
%PE 9/23/2022 - Cleaned for final release
%PE 5/24/2022 - Added additional parameters for PV_ISO calculations
%PE 4/28/2022 - Adapted for multiple methods/fields (center and contour)
%PE 4/27/2022 - Adapted from original extraction_eddies_new
%--------------------------------------
%Purpose: Grabs local maxes/mins and contours around those maxes/mins, classifies those as eddies
%Inputs: X (vector of lons), Y (vector of lats), centerField (field to be used for local extrema)
%        edgeField (field to be used for contours), minDist (min distance between eddies),
%        minAmp (min amp of eddy center), minRad (min radius of eddy), minRatio (minimum shape
%        parameter), EKE (EKE), centerMethod ("SSH", "PV_ISO", "edgeField{1,5}", "MV", "OW", "LNAM"), edgeMethod
%        ("SSH", "PV_ISO", "edgeField{1,5}", "MV", "OW", "LNAM")
%Outputs: Subsurface Eddies of the style chosen
%--------------------------------------
function [Anticyclonic_Cell,Cyclonic_Cell,Nanti,Ncyclo,Label_anti,Label_cyclo,maxRet] = Paul_OPT_EddyExtractionMaster(X,Y,centerField,edgeField,minDist,minAmp,minRad,minRatio,EKE,centerMethod,edgeMethod)
%% Configure inputs
% Create matrices of coordinate

X = double(X); Y = double(Y);
[ygrid,xgrid] = meshgrid(Y,X);

% Initialize number of eddies at 0
Ncyclo = 0;
Nanti = 0;

% Critical variable for window size in lat/lon (5 degrees is sufficient for all mesoscale eddies)
maxWindow = 5;

% SET SMOOTHING IF NEEDED
if ((centerMethod == "VORT_WA") || (centerMethod == "OW_WA") || (centerMethod == "PV_ISO"))
    if centerField{1,2} == 0
        centerField{1,1} = centerField{1,1};
    else
        centerField{1,1} = spSmoothing(centerField{1,1},centerField{1,2});
    end
    centerField = centerField{1,1};
end
if ((centerMethod == "VORT_T") || (centerMethod == "OW_T"))
    if centerField{1,4} == 0
        centerField{1,1} = centerField{1,1};
    else
        centerField{1,1} = spSmoothing(centerField{1,1},centerField{1,4});
    end
end

%% Find centers
% Look for non NAN local maximum (Anticylone) and minimum (Cyclone)
% Store them in a sorted vectors maxMat without NaN

% WINDING ANGLE ALGORITHMS
if ((centerMethod == "SSH") || (centerMethod == "PV_ISO") || (centerMethod == "VORT_WA")  || (centerMethod == "OW_WA"))

    f1 = localMaximum(centerField,minDist);  f1(isnan(centerField(f1)))=[];
    f2 = localMaximum(-centerField,minDist);  f2(isnan(centerField(f2)))=[];
    %Removing positives
    if (centerMethod == "OW_WA")
        f2 = [];
    end
    maxMat=sort([f1(:);f2(:)]);
    maxRet={f1,f2};

    % THRESHOLD ALGORITHMS
elseif ((centerMethod == "VORT_T")  || (centerMethod == "OW_T"))
    %Removing positives
    if (centerMethod == "OW_T")
        centerField{1,1}(centerField{1,1}<0) = 0;
    end
    [maxMat, maxRet] = Threshold_eddy_centers(centerField{1,1},centerField{1,2},X,Y,centerField{1,3});
    centerField = centerField{1,1};

    % AMEDA ALGORITHM
elseif ((centerMethod == "LNAM") || (centerMethod == "LNAM_ISO"))

    [maxMat, maxRet] = AMEDA_eddy_centers(centerField{1,1},centerField{1,2},X,Y,centerField{1,3},centerField{1,4},centerField{1,5});
    centerField = centerField{1,1};

end

%% Create Cells for Cyclone and Anticylone. Row dimension large but will be
% reduced when cells will be filled
Cyclonic_Cell=cell(length(maxMat),12);
Anticyclonic_Cell=cell(length(maxMat),12);

%% Time step on each center ID'd earlier - Winding Angle (i.e. closed contour of edgeField)
if ((edgeMethod == "SSH") || (edgeMethod == "VORT_WA") || (edgeMethod == "PV_ISO") || (edgeMethod == "OW_WA"))
    %If we need to set smoothing
    if ((edgeMethod == "VORT_WA") || (edgeMethod == "OW_WA") || (edgeMethod == "PV_ISO"))
        if edgeMethod == "OW_WA"
            edgeField{1,1}(edgeField{1,1}<0) = NaN; %remove positive values
            edgeField{1,3}(edgeField{1,3}>0) = 1; edgeField{1,3}(edgeField{1,3}<0) = -1; %get vort scaling
            edgeField{1,1} = abs(edgeField{1,1}).*edgeField{1,3}; %rescale for AE/CE
        end
        if edgeField{1,2} == 0
            edgeField = edgeField{1,1};
        else
            edgeField = spSmoothing(edgeField{1,1},edgeField{1,2});
        end
    end
    %The "thresh" value for contour searching - not sensitive and typically lower than minimum
    switch edgeMethod
        case "SSH"
            thresh = 1e-4;
        case "VORT_WA"
            thresh = 1e-6;
        case "OW_WA"
            thresh = 1e-14;
        case "PV_ISO"
            thresh = 1e-7;
    end

    %Loop on local extreme
    for i=1:length(maxMat)

        % thisMax is the local SSH maximum considered
        thisMax = edgeField(maxMat(i));

        %Ensure compatibility
        if isnan(thisMax)
            continue
        end

        % Find indices of longitudes and lattitude around maximum +/- maxWindow
        % degrees. For vectors (needed for contours function)
        % and in the matrices
        indlon = find(X>=xgrid(maxMat(i))-maxWindow & X<=xgrid(maxMat(i))+maxWindow);
        indlat = find(Y<=ygrid(maxMat(i))+maxWindow & Y>=ygrid(maxMat(i))-maxWindow);
        indf =  find(ygrid(maxMat)<=ygrid(maxMat(i))+maxWindow &...
            ygrid(maxMat)>=ygrid(maxMat(i))-maxWindow &...
            xgrid(maxMat)>=xgrid(maxMat(i))-maxWindow &...
            xgrid(maxMat)<=xgrid(maxMat(i))+maxWindow);

        %Generate local grid
        localField = edgeField(indlon,indlat);
        xxlocal = xgrid(indlon,indlat);
        yylocal = ygrid(indlon,indlat);

        % Binary search algorithm (Dichotomie) to find contour of eddies [PEGLIASCO CODE]

        % Anticyclones --------------------------------------------------------------------------------

        % Set initial amplitude differences with center of 2 m (fin-deb)
        % finish2=fin at beginnig. IsGood is defined=0 as no contour has be
        % identified yet.
        deb = 0;
        fin = 1+thresh;
        finish2 = fin;
        IsGood = 0;

        % iind will be used to store the anplitude and index of identified
        % contour corresponding to criterion in the area studied
        iind = [];

        % While loop on contours espaced by 1 mm
        while fin-deb>thresh

            % Store in c all contour with values hd-fin in the perimeter
            % defined presviously indicated by indlat && indlon
            c=contourc(X(indlon),Y(indlat),edgeField(indlon,indlat)',[thisMax-fin thisMax-fin]);

            % If there is at least one contour detected
            if isempty(c)==0
                % Split contours in a structure array of contours
                [cstruct c] = splitcontours(c);
                % Count the number of contours stored in the structure
                NbreC = length(cstruct);
                % Loop on these contours
                for j=1:NbreC
                    % Store contours X and Y postion in vectors xd and yd
                    xd = cstruct(j).x;
                    yd = cstruct(j).y;
                    %plot(xd,yd,'r-')
                    % If it is a closed contour
                    if xd(end) == xd(1) && yd(end) ==yd(1)
                        %Testing ratios
                        geom = polygeom(xd,yd);
                        centroidX = geom(2); centerX = xgrid(maxMat(i));
                        centroidY = geom(3); centerY = ygrid(maxMat(i));
                        distMat = zeros(2,length(xd));
                        for lenCont = 1:length(xd)
                            distMat(1,lenCont) = sqrt((xd(lenCont)-centroidX)^2 + (yd(lenCont)-centroidY)^2);
                            distMat(2,lenCont) = sqrt((xd(lenCont)-centerX)^2 + (yd(lenCont)-centerY)^2);
                        end
                        minorAxisOID = min(distMat(1,:));
                        majorAxisOID = max(distMat(1,:));
                        ratioAxesOID = minorAxisOID/majorAxisOID;
                        minorAxis = min(distMat(2,:));
                        majorAxis = max(distMat(2,:));
                        ratioAxes = minorAxis/majorAxis;
                        % If the local maximum lied in this contour and if
                        % there is only on loca maximum using mex inpoly.c
                        %figure; plot(xd,yd); hold on; scatter(xgrid(maxMat),ygrid(maxMat));
                        if inpoly([xgrid(maxMat(i));ygrid(maxMat(i))],[xd(:) yd(:)]')==1 && ...
                                length(find(inpoly([xgrid(maxMat(indf))';ygrid(maxMat(indf))'],[xd(:) yd(:)]'))==1)==1 &&...
                                ((ratioAxesOID>minRatio) && (ratioAxes>minRatio/2)) %ratiotest
                            iiind = find(inpoly([xxlocal(:) yylocal(:)]',[xd(:) yd(:)]')==1);
                            if all(isfinite(localField(iiind)))==1 %& length(isfinite(ssh(iiind)))>=4
                                % Give flag of 1 to the variable is good that
                                % saying that a contour corresponds to defined
                                % criterions
                                IsGood = 1;
                                iind = [fin j];
                                break
                            end
                        end
                    end
                end
            end

            % Condition on the detection of contour corresponding to criterions
            if IsGood == 1 % One contour was detected

                % Divide and shift the area of research toward exterior
                % (toward fin). ie Divide distance between contour and
                % exterior by 2. Reset IsGood to 0.
                deb = fin;
                fin = deb + (finish2-deb)/2;
                IsGood = 0;

            elseif IsGood == 0 % No contour was detected

                % Divide and shift the area of research toward interior
                % (toward center). ie Divide distance between contour and
                % center by 2.
                finish2 = fin;
                fin = deb + (finish2-deb)/2;
            end
        end

        % If a contour stored in ind was detected with amplitude superior
        % than the minimum defined when the function was called
        if isempty(iind)==0 & iind(1)>=minAmp &...
                length(isfinite(localField(iiind)))>=4 &...
                length(unique(xxlocal(iiind)))>=2 &...
                length(unique(yylocal(iiind)))>=2

            % Find contours delimiting the closed area around one center
            % where the difference of amplitude with the local maximum is
            % maximum and that correspond to previous criterions
            c=contourc(X(indlon),Y(indlat),edgeField(indlon,indlat)',[thisMax-iind(1) thisMax-iind(1)]);

            % Split contours and take the contours values identified by the
            % iind(2) which store index of good contour (j)
            [cstruct ~] = splitcontours(c);
            xd = cstruct(iind(2)).x;
            yd = cstruct(iind(2)).y;

            [b,a] = calcul_aire_rayon(xd,yd);

            % Calculate the ratio parameters (distance from center, distance from centroid)
            geom = polygeom(xd,yd);
            centroidX = geom(2); centerX = xgrid(maxMat(i));
            centroidY = geom(3); centerY = ygrid(maxMat(i));
            distMat = zeros(2,length(xd));
            for lenCont = 1:length(xd)
                distMat(1,lenCont) = sqrt((xd(lenCont)-centroidX)^2 + (yd(lenCont)-centroidY)^2);
                distMat(2,lenCont) = sqrt((xd(lenCont)-centerX)^2 + (yd(lenCont)-centerY)^2);
            end
            minorAxisOID = min(distMat(1,:));
            majorAxisOID = max(distMat(1,:));
            ratioAxesOID = minorAxisOID/majorAxisOID;
            minorAxis = min(distMat(2,:));
            majorAxis = max(distMat(2,:));
            ratioAxes = minorAxis/majorAxis;

            %Check radius and ratio, both for center and centroid (center is less stringent)
            if ((a>=minRad) && ((ratioAxesOID>minRatio) && (ratioAxes>minRatio/2)))

                % Indicates that one more anticylones was detected
                Nanti = Nanti+1;

                % Fill values of X and Y center
                Anticyclonic_Cell{Nanti,1}=xgrid(maxMat(i));
                Anticyclonic_Cell{Nanti,2}=ygrid(maxMat(i));

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
                Anticyclonic_Cell{Nanti,9}=single(iind(1));

                % Fill ratio values
                Anticyclonic_Cell{Nanti,11}=single(ratioAxesOID);
                Anticyclonic_Cell{Nanti,12}=single(ratioAxes);

                %TEST
                %             else
                %                 figure; plot(xd, yd); hold on;
                %                 scatter(xgrid(maxMat(i)), ygrid(maxMat(i)), 'b');
                %                 scatter(single(geom(2)), single(geom(3)), 'r');
            end

        end

        % Cyclones ------------------------------------------------------------------------------------

        % Set initial amplitude differences with center of 200 cm (fin-deb)
        % finish2=fin at beginnig. IsGood is defined=0 as no contour has be
        % identified yet.
        deb = 0;
        fin = 1+thresh;
        finish2 = fin;
        IsGood = 0;

        % iind will be used to store the anplitude and index of identified
        % contour corresponding to criterion in the area studied
        iind = [];

        % While loop on contours espaced by 1 mm
        while fin-deb>thresh

            % Store in c all contour with values hd-fin in the perimeter
            % defined presviously indicated by indlat && indlon
            c=contourc(double(X(indlon)),double(Y(indlat)),edgeField(indlon,indlat)',[thisMax+fin thisMax+fin]);

            % If there is at least one contour detected
            if isempty(c)==0
                % Split contours in a structure array of contours
                [cstruct ~] = splitcontours(c);
                % Count the number of contours stored in the structure
                NbreC = length(cstruct);
                % Loop on these contours
                for j=1:NbreC
                    % Store contours X and Y postion in vectors xd and yd
                    xd = cstruct(j).x;
                    yd = cstruct(j).y;
                    % If it is a closed contour
                    if xd(end) == xd(1) && yd(end) ==yd(1)
                        %Testing ratios
                        geom = polygeom(xd,yd);
                        centroidX = geom(2); centerX = xgrid(maxMat(i));
                        centroidY = geom(3); centerY = ygrid(maxMat(i));
                        distMat = zeros(2,length(xd));
                        for lenCont = 1:length(xd)
                            distMat(1,lenCont) = sqrt((xd(lenCont)-centroidX)^2 + (yd(lenCont)-centroidY)^2);
                            distMat(2,lenCont) = sqrt((xd(lenCont)-centerX)^2 + (yd(lenCont)-centerY)^2);
                        end
                        minorAxisOID = min(distMat(1,:));
                        majorAxisOID = max(distMat(1,:));
                        ratioAxesOID = minorAxisOID/majorAxisOID;
                        minorAxis = min(distMat(2,:));
                        majorAxis = max(distMat(2,:));
                        ratioAxes = minorAxis/majorAxis;
                        % If the local maximum lied in this contour and if
                        % there is only on loca maximum (Using mex
                        if inpoly([xgrid(maxMat(i));ygrid(maxMat(i))],[xd(:) yd(:)]')==1 &&...
                                length(find(inpoly([xgrid(maxMat(indf))';ygrid(maxMat(indf))'],[xd(:) yd(:)]'))==1)==1 &&...
                                ((ratioAxesOID>minRatio) && (ratioAxes>minRatio/2)) %ratiotest
                            iiind = find(inpoly([xxlocal(:) yylocal(:)]',[xd(:) yd(:)]')==1);
                            if all(isfinite(localField(iiind)))==1

                                % Give flag of 1 to the variable is good that
                                % saying that a contour corresponds to defined
                                % criterions
                                IsGood = 1;
                                iind = [fin j];
                                break
                            end
                        end
                    end
                end
            end

            % Condition on the detection of contour corresponding to criterions
            if IsGood == 1 % One contour was detected

                % Divide and shift the area of research toward exterior
                % (toward fin). ie Divide distance between contour and
                % exterior by 2. Reset IsGood to 0.
                deb = fin;
                fin = deb + (finish2-deb)/2;
                IsGood = 0;

            elseif IsGood == 0 % No contour was detected

                % Divide and shift the area of research toward interior
                % (toward center). ie Divide distance between contour and
                % center by 2.
                finish2 = fin;
                fin = deb + (finish2-deb)/2;
            end
        end

        % If a contour stored in ind was detected with amplitude superior
        % than the minimum defined when the function was called
        if isempty(iind)==0 & iind(1)>=minAmp &...
                length(isfinite(localField(iiind)))>=4 &...
                length(unique(xxlocal(iiind)))>=2 &...
                length(unique(yylocal(iiind)))>=2

            % Find contours delimiting the closed area around one center
            % where the difference of amplitude with the local maximum is
            % maximum and that correspond to previous criterions
            c=contourc(double(X(indlon)),double(Y(indlat)),edgeField(indlon,indlat)',[thisMax+iind(1) thisMax+iind(1)]);

            % Split contours and take the contours values identified by the
            % iind(2) which store index of good contour (j)
            [cstruct ~] = splitcontours(c);
            xd = cstruct(iind(2)).x;
            yd = cstruct(iind(2)).y;

            [b,a] = calcul_aire_rayon(xd,yd);

            % Calculate the ratio parameters (distance from center, distance from centroid)
            geom = polygeom(xd,yd);
            centroidX = geom(2); centerX = xgrid(maxMat(i));
            centroidY = geom(3); centerY = ygrid(maxMat(i));
            distMat = zeros(2,length(xd));
            for lenCont = 1:length(xd)
                distMat(1,lenCont) = sqrt((xd(lenCont)-centroidX)^2 + (yd(lenCont)-centroidY)^2);
                distMat(2,lenCont) = sqrt((xd(lenCont)-centerX)^2 + (yd(lenCont)-centerY)^2);
            end
            minorAxisOID = min(distMat(1,:));
            majorAxisOID = max(distMat(1,:));
            ratioAxesOID = minorAxisOID/majorAxisOID;
            minorAxis = min(distMat(2,:));
            majorAxis = max(distMat(2,:));
            ratioAxes = minorAxis/majorAxis;

            %Check radius and ratio, both for center and centroid
            if ((a>=minRad) && ((ratioAxesOID>minRatio) && (ratioAxes>minRatio/2)))
                % Indicates that one more anticylones was detected
                Ncyclo = Ncyclo+1;

                % Fill values of X and Y center
                Cyclonic_Cell{Ncyclo,1}=xgrid(maxMat(i));
                Cyclonic_Cell{Ncyclo,2}=ygrid(maxMat(i));

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
                Cyclonic_Cell{Ncyclo,9}=single(iind(1));

                %Fill ratios
                Cyclonic_Cell{Ncyclo,11}=single(ratioAxesOID);
                Cyclonic_Cell{Ncyclo,12}=single(ratioAxes);

                %DEBUG PLOT
                %else
                %figure; plot(xd, yd); hold on;
                %scatter(xgrid(maxMat(i)), ygrid(maxMat(i)), 'b');
                %scatter(single(geom(2)), single(geom(3)), 'r');

            end

        end

        %Reduce cells to their non-NaN extend
        Anticyclonic_Cell(Nanti+1:end,:) = [];
        Cyclonic_Cell(Ncyclo+1:end,:) = [];
    end

    %% Time step on each center ID'd earlier: THRESHOLD METHODS
elseif ((edgeMethod == "VORT_T") || (edgeMethod == "OW_T"))
    %Smoothing
    if edgeField{1,4} == 0
        edgeField{1,1} = edgeField{1,1};
    else
        edgeField{1,1} = spSmoothing(edgeField{1,1},edgeField{1,4});
    end
    if edgeMethod == "OW_T"
        edgeField{1,1}(edgeField{1,1}<0) = NaN; %remove positive values
        edgeField{1,5}(edgeField{1,5}>0) = 1; edgeField{1,5}(edgeField{1,5}<0) = -1; %get vort scaling
        edgeField{1,1} = abs(edgeField{1,1}).*edgeField{1,5}; %rescale for AE/CE
    end
    [Anticyclonic_Cell, Cyclonic_Cell, Ncyclo, Nanti] = Threshold_eddy_edges(edgeField{1,1},edgeField{1,2},X,Y,edgeField{1,3},maxMat,minRad,minAmp,minRatio);
    edgeField = edgeField{1,1};

end

%% Reverse Cyclonic and Anticyclonic values if edgeField{1,5}/PV_ISO (because is inverse to SSH signs of cyclone anticyclone, WATCH HEMISPHERE)
if ((edgeMethod == "VORT_WA") || (edgeMethod == "PV_ISO") || (edgeMethod == "VORT_T") || (edgeMethod == "OW_WA") || (edgeMethod == "OW_T"))
    Intermed = Anticyclonic_Cell;
    Anticyclonic_Cell = Cyclonic_Cell;
    Cyclonic_Cell = Intermed;
    Intermed = Nanti;
    Nanti = Ncyclo;
    Ncyclo = Intermed;
end

%% LABEL MATRICES
Label_anti = single(zeros(length(X),length(Y)));
Label_cyclo = single(zeros(length(X),length(Y)));
[y,x] = meshgrid(Y,X);
if Ncyclo>0
    for j=1:Ncyclo %Loop on each Cyclone
        % Edge coordinates
        disp(['CE loop:   ' int2str(j)  ' of '  int2str(Ncyclo)]);
        Xcon = double(Cyclonic_Cell{j,5});
        Ycon = double(Cyclonic_Cell{j,6});

        % Change x && y to reduce search area
        ind1 = find(x>=min(Xcon)-0.5 & x<=max(Xcon)+0.5 & y>=min(Ycon)-0.5 & y<=max(Ycon)+0.5);

        % Coordinate index of grid points inside the eddy
        in = find(inpoly([x(ind1),y(ind1)]',[Xcon(:),Ycon(:)]')==1);
        INDinside = ind1(in);

        % Xcenter and Ycenter
        [~,iiind] = min(centerField(INDinside));
        INDcenter = INDinside(iiind(1));

        Cyclonic_Cell{j,1}= single(x(INDcenter));
        Cyclonic_Cell{j,2}= single(y(INDcenter));

        % Mean EKE
        Cyclonic_Cell{j,10}=single(nanmean(EKE(INDinside)));

        Label_cyclo(ind1(in)) = j;
    end
end

% Compute and store parameters of each Anticyclones if there is at least
% one detected
if Nanti>0
    for j=1:Nanti %Loop on each Cyclone
        disp(['AE loop:   ' int2str(j)  ' of '  int2str(Nanti)]);
        % Edge coordinates
        Xcon = double(Anticyclonic_Cell{j,5});
        Ycon = double(Anticyclonic_Cell{j,6});

        % Change x && y to reduce search area
        ind1 = find(x>=min(Xcon)-0.5 & x<=max(Xcon)+0.5 & y>=min(Ycon)-0.5 & y<=max(Ycon)+0.5);

        % Coordinate index of grid points inside the eddy
        a1=[x(ind1),y(ind1)]';b1=[Xcon(:),Ycon(:)]';%pause(20);
        in = find(inpoly(a1,b1)==1);%pause(5)
        INDinside = ind1(in);

        % Xcenter and Ycenter
        [~,iiind] = max(centerField(INDinside));
        INDcenter = INDinside(iiind(1));

        Anticyclonic_Cell{j,1}= single(x(INDcenter));
        Anticyclonic_Cell{j,2}= single(y(INDcenter));

        % Mean EKE
        Anticyclonic_Cell{j,10}=single(nanmean(EKE(INDinside)));

        Label_anti(ind1(in)) = j;
    end
end




