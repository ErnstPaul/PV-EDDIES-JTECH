%% Paul_OPT_EdgeSmoothingAndRatio.m (version 3.1)
%Author: Paul Ernst
%Date Created: 5/10/2022
%Date of Last Update: 1/12/2023
%Update History:
%PE 1/12/2023 - Cleared in final revision
%PE 6/2/2022 - Adjusted to work for EdgeSmoothing
%PE 5/24/2022 - Adjusted for PV_ISO corrections from Yves Morel
%PE 5/11/2022 - Added SSH detection and altered params
%PE 5/10/2022 - Created
%--------------------------------------
%Purpose: Takes a request for comparing eddy tracking methods and does so (EDGESMTHRAT, FIG4)
%Inputs: All defaults, individual files created by file preparation steps
%        The contour/edge methods we are trying are:
%           0: SSH (winding angle, ground truth method) (CALLED SSH)
%           1: Isopycnal PV (CALLED PV_ISO, smoothed by factor smoothingFactor)
%           2: Vorticity (CALLED VORT_T and VORT_WA, smoothed by factor smoothingFactor with stdFactor)
%           3: OW (raw) (CALLED OW_T and OW_WA, smoothed by factor smoothingFactor with stdFactor)
%Outputs: Comparison table and separate files for each algorithm
%--------------------------------------
%% Inputs
function Paul_OPT_EdgeSmoothingAndRatio(inputDir)
load('DefaultData_OPT.mat') %See Paul_OPT_CreateDefaultSettings.m
Paul_OPT_LoadDefaultSettings(); %See corresponding function

%Title for info mat
titlestrings = ["OPT_E_VORT_T_SMTHRAT" "OPT_E_VORT_WA_SMTHRAT" "OPT_E_OW_T_SMTHRAT" "OPT_E_OW_WA_SMTHRAT" "OPT_E_PV_ISO_SMTHRAT"];
RPT = length(titlestrings);
%METHODS
centerMethods = ["SSH"]; %list of methods to find center: SSH, PV_ISO, VORT, MV, OW, VG, LNAM
edgeMethodss = ["VORT_T","VORT_WA","OW_T","OW_WA","PV_ISO"]; %list of methods to find edge: SSH, PV_ISO, VORT, MV, OW, VG

%ONLY USE FOR OPTIMIZING BOTH SMOOTH AND RATIO WITH MANY METHODS %CHANGE THIS BACK TO 1
for RPTLOOP = 1:RPT
    if (RPT > 1)
        titlestring = titlestrings(RPTLOOP);
        edgeMethods = edgeMethodss(RPTLOOP);
    end
    %ONLY USE FOR OPTIMIZING BOTH SMOOTH AND RATIO WITH MANY METHODS
    createPairs = true; %If true, then creates a total matchup of all pairings; if false, takes only pairs above
    K = 0.95; stdFactor = 3;

    %GRAB PREVIOUS DATA FOR OPTIMIZATION
    if (edgeMethods == "VORT_T")
        prevString = "OPT_E_VORT_T_STD";
        load("IO_" + prevString + "_Comparisons.mat", "winner")
        stdFactor = winner;
    elseif (edgeMethods == "OW_T")
        prevString = "OPT_E_OW_T_STD";
        load("IO_" + prevString + "_Comparisons.mat", "winner")
        stdFactor = winner;
    end

    %Factor for thresholds
    smoothingFactor = 0:15; %0 - 15
    a = 10; b = 8; %0 -10, factors for VG
    minAmp = 0;
    minRatio = 0:0.02:0.3; %0 - 0.3 is our minimum threshold for ellipsoidness
    minRad = 25; %minimum radius (we're looking at mesoscale only, so 25+)
    minDist = 25;

    %Logicals for if we're optimizing something, pick only one of the following:
    optimizeSTD = false;
    optimizeK = false;
    optimizeSmooth = false;
    optimizeMethod = false;
    optimizeDist = false;
    optimizeRatio = false;
    optimizeSmoothAndRatio = true;

    %INPUT FILES
    % - SEE Paul_OPT_Step1_ComparisonsMO_SS.m - First name of "pv"
    %These files will have:
    %Isopycnally Averaged PV_ISO at Surface
    %Isopycnally Averaged VORT at Surface
    %Surface: U, V, S, T, RHO, O-W Parameter, SSH, Vorticity, LNAM, EKE, X and Y

    %OUTPUT FILES
    %Location we're sending each of our file lists out to
    outputDir = strcat(basepath, '/COMPARISON/');
    [status,message,messageid] = mkdir(outputDir);

    %KEY FILES
    %Location of winding angle files to compare results to (searches by date)
    keyDir = strcat(basepath, '/COMPARISON/KEYFILES/');

    %% Enable file reading
    %Filter for usable files in the input directory
    list = dir([ inputDir '/' 'pv*.mat']);
    configFile = dir([ inputDir '/' 'Config*.mat']);
    load([configFile.folder '/' configFile.name]);

    %% Construct master matrix of inputs
    %If true, matches up each method; if false, only matches what we have as an input
    if optimizeMethod
        if createPairs
            masterMethods = strings(2,length(centerMethods)*length(edgeMethods));
            count = 1;
            for i = 1:length(centerMethods)
                for j = 1:length(edgeMethods)
                    %This gets every pair as (A | B) =/= (B | A)
                    masterMethods(1,count) = centerMethods(1,i);
                    masterMethods(2,count) = edgeMethods(1,j);
                    count = count + 1;
                end
            end
        else
            masterMethods = [centerMethods; edgeMethods];
        end
        [~, optLen] = size(masterMethods);
        masterDirList = strings(length(list),optLen);
    else
        masterMethods = [centerMethods; edgeMethods];
    end

    %% Sub-settings prior to file loop
    %Finding minimum radius/window size based on resolutions
    filename = [inputDir '/' list(1).name];
    load(filename,'X','Y')
    Rmin = ceil(2*max([diff(unique(X(:)));diff(unique(Y(:)))])*111.1);
    disp(['Considering the actual resolution and in order to include']);
    disp(['at least 4 grid points inside the detected eddies,'])
    disp(['the minimum eddy radius will be of ' num2str(Rmin) ' km']);
    disp('')
    minRad = max([minRad Rmin]);
    disp('')
    disp(['The windows size is ' , num2str(minDist*(X(2)-X(1))), ' degrees, corresponding to ',num2str(minDist), ' pixels: '])

    %% Loop on each file

    %determine non-method optimization
    if optimizeK
        optLen = length(K);
    end
    if optimizeSmooth
        optLen = length(smoothingFactor);
    end
    if optimizeSTD
        optLen = length(stdFactor);
    end
    if optimizeDist
        optLen = length(minDist);
    end
    if optimizeRatio
        optLen = length(minRatio);
    end
    if optimizeSmoothAndRatio
        count = 1;
        masterSR = zeros(2,length(smoothingFactor)*length(minRatio));
        optLen = length(smoothingFactor)*length(minRatio);
        for i = 1:length(smoothingFactor)
            for j = 1:length(minRatio)
                %This gets every pair as (A | B) =/= (B | A)
                masterSR(1,count) = smoothingFactor(1,i);
                masterSR(2,count) = minRatio(1,j);
                count = count + 1;
            end
        end
    end

    methCount = 1;

    totalSims = zeros(1,optLen);
    totalRatio = zeros(1,optLen);
    totalSpatialR = zeros(1,optLen);
    totalSpatialC = zeros(1,optLen);
    totalNum = zeros(1,optLen);
    %Configure center method

    optStart = 1;

    for optN = optStart:optLen
        if optimizeMethod
            centerMethod = masterMethods(1,optN);
            %Configure edge method
            edgeMethod = masterMethods(2,optN);
        else
            centerMethod = centerMethods;
            edgeMethod = edgeMethods;
        end

        %Configure output directory and create if needed
        thisOutputDir = [outputDir + "/" + centerMethod + "_" + edgeMethod];

        %Adjust parameters based on optimization
        if optimizeK
            thisK = K(optN);
            thisOutputDir = [thisOutputDir + optN + "K"];
            disp(["Optimizing K: " + num2str(thisK)])
        else
            thisK = K;
        end
        if optimizeSmooth
            thisSmoothingFactor = smoothingFactor(optN);
            thisOutputDir = [thisOutputDir + optN + "SM"];
            disp(["Optimizing Smoothing Factor: " + num2str(thisSmoothingFactor)])
        end
        if optimizeSTD
            thisSTD = stdFactor(optN);
            thisOutputDir = [thisOutputDir + optN + "STD"];
            disp(["Optimizing Standard Deviation Modifier: " + num2str(thisSTD)])
        else
            thisSTD = stdFactor;
        end
        if optimizeDist
            thisDist = minDist(optN);
            thisOutputDir = [thisOutputDir + optN + "DST"];
            disp(["Optimizing Minimum Distance: " + num2str(thisDist)])
        else
            thisDist = minDist;
        end
        if optimizeRatio
            thisRatio = minRatio(optN);
            thisOutputDir = [thisOutputDir + optN + "RAT"];
            disp(["Optimizing Minimum Ratio: " + num2str(thisRatio)])
        end
        if optimizeSmoothAndRatio
            thisSmoothingFactor = masterSR(1,optN);
            thisRatio = masterSR(2,optN);
            disp(["Optimizing Minimum Ratio: " + num2str(thisRatio)])
            disp(["Optimizing Smoothing Factor: " + num2str(thisSmoothingFactor)])
            thisOutputDir = [thisOutputDir + optN + "SMTHRAT"];
        elseif ((~optimizeRatio) && (~optimizeSmooth))
            thisSmoothingFactor = smoothingFactor;
            thisRatio = minRatio;
        elseif (~optimizeRatio)
            thisRatio = minRatio;
        elseif (~optimizeSmooth)
            thisSmoothingFactor = smoothingFactor;
        end

        %If this is all SSH
        if (strcmp(centerMethod,"SSH") && strcmp(edgeMethod,"SSH"))
            thisOutputDir = [outputDir + "/" + "KEYFILES/"];
        end
        [status,message,messageid] = mkdir(thisOutputDir);
        masterDirList(methCount) = thisOutputDir;

        disp("Beginning eddy extraction: " + centerMethod + " Centers with " + edgeMethod + " Edges.")

        %Configure list
        for fileNum = 1:length(list)
            tic
            disp(['Extracting eddies for:    ' list(fileNum).name ''])
            filename = [inputDir  list(fileNum).name];
            %Load and reaffirm presence of all necessary variables
            load(filename)
            EKE = double(EKE(:,:,1));
            U = double(U(:,:,1));
            V = double(V(:,:,1));
            PV_ISO = double(PV_ISO(:,:));
            VORT = double(VORT(:,:,1));
            OW = double(OW(:,:));
            LOW = double(LOW(:,:,1));
            ZOS = double(ZOS(:,:));
            LNAM = double(LNAM(:,:,1));

            %Center detection switch cases
            switch centerMethod
                case "SSH"
                    centerField = ZOS;
                case "PV_ISO"
                    centerField = cell(1,2);
                    centerField{1,1} = PV_ISO;
                    centerField{1,2} = thisSmoothingFactor;
                case "VORT_WA"
                    centerField = cell(1,2);
                    centerField{1,1} = VORT;
                    centerField{1,2} = thisSmoothingFactor;
                case "VORT_T"
                    centerField = cell(1,4);
                    centerField{1,1} = VORT;
                    centerField{1,2} = VORTSTATS;
                    centerField{1,3} = thisSTD;
                    centerField{1,4} = thisSmoothingFactor;
                case "OW_WA"
                    centerField = cell(1,2);
                    centerField{1,1} = OW;
                    centerField{1,2} = thisSmoothingFactor;
                case "OW_T"
                    centerField = cell(1,4);
                    centerField{1,1} = OW;
                    centerField{1,2} = OWSTATS;
                    centerField{1,3} = thisSTD;
                    centerField{1,4} = thisSmoothingFactor;
                case "LNAM"
                    centerField = cell(1,5);
                    centerField{1,1} = LOW;
                    centerField{1,2} = LNAM;
                    centerField{1,3} = U;
                    centerField{1,4} = V;
                    centerField{1,5} = thisK;
            end

            %Edge detection switch cases
            switch edgeMethod
                case "SSH"
                    edgeField = ZOS;
                case "PV_ISO"
                    edgeField = cell(1,2);
                    edgeField{1,1} = PV_ISO;
                    edgeField{1,2} = thisSmoothingFactor;
                case "VORT_WA"
                    edgeField = cell(1,2);
                    edgeField{1,1} = VORT;
                    edgeField{1,2} = thisSmoothingFactor;
                case "VORT_T"
                    edgeField = cell(1,4);
                    edgeField{1,1} = VORT;
                    edgeField{1,2} = VORTSTATS;
                    edgeField{1,3} = thisSTD; %thisSTD
                    edgeField{1,4} = thisSmoothingFactor;
                case "OW_WA"
                    edgeField = cell(1,3);
                    edgeField{1,1} = OW;
                    edgeField{1,2} = thisSmoothingFactor;
                    edgeField{1,3} = VORT; %to distinguish AE (neg) from CE (pos)
                case "OW_T"
                    edgeField = cell(1,5);
                    edgeField{1,1} = OW;
                    edgeField{1,2} = OWSTATS;
                    edgeField{1,3} = thisSTD; %thisSTD
                    edgeField{1,4} = thisSmoothingFactor;
                    edgeField{1,5} = VORT; %to distinguish AE (neg) from CE (pos)
            end

            %Run algorithm that we've constructed
            [Anticyclonic_Cell,Cyclonic_Cell,Nanti,Ncyclo,Label_anti,Label_cyclo,centerLabels] = ...
                Paul_OPT_EddyExtractionMaster(X,Y,centerField,edgeField,thisDist,...
                minAmp,minRad,thisRatio,EKE,centerMethod,edgeMethod);
            FieldsCells = [
                '1.Longitude of the Eddy Center          '; ...
                '2.Latitude of the Eddy Center           '; ...
                '3.Longitude of the Eddy Centroid        '; ...
                '4.Latitude of the Eddy Centroid         '; ...
                '5.Longitudes of the Eddy Edge           '; ...
                '6.Latitudes of the Eddy Edge            '; ...
                '7.Equivalent Radius [km]                '; ...
                '8.Eddy Area [km^2]                      '; ...
                '9.Amplitude [m]                         '; ...
                '10.Mean EKE [(m/s)^2]                   '; ...
                '11.Ratio of Centroid Axes               '; ...
                '12.Ratio of Center Axes                 '];

            %Save initial matrix back into output, adaptive based on where it is
            date2use = datestr(date_num,'mm/dd/yyyy');
            date2use = [date2use(7:end) date2use(1:2) date2use(4:5)];
            filename_out = [thisOutputDir + "/comp_" + date2use];
            save(filename_out,'Anticyclonic_Cell','Cyclonic_Cell','Nanti','Ncyclo','centerLabels',...
                'Label_anti','Label_cyclo','ZOS','centerField','edgeField','date_num','FieldsCells');

            %No similarity to be done if everything is SSH
            if (strcmp(centerMethod,"SSH") && strcmp(edgeMethod,"SSH"))
                continue;
            end

            %Run similarity calculation
            %Grab the date
            filteredFile = dir([ keyDir '/' 'comp_' date2use '.mat']);
            keyname = [keyDir  filteredFile(1).name];
            %Before loading, rename our first variables
            Anticyclonic_Cell2 = Anticyclonic_Cell;
            Cyclonic_Cell2 = Cyclonic_Cell;
            Nanti2 = Nanti;
            Ncyclo2 = Ncyclo;
            Label_anti2 = Label_anti;
            Label_cyclo2 = Label_cyclo;
            %Now Load
            load(keyname)

            %Get logicals
            Label_anti1 = double(Label_anti>0); Label_cyclo1 = double((Label_cyclo>0)); %key
            Label_anti2 = double(Label_anti2>0); Label_cyclo2 = double(Label_cyclo2>0); %test

            %Get ratio matrices
            ratioKey = nanmean(cell2mat([Anticyclonic_Cell(:,11); Anticyclonic_Cell(:,12); Cyclonic_Cell(:,11); Cyclonic_Cell(:,12)]));
            ratioTest = nanmean(cell2mat([Anticyclonic_Cell2(:,11); Anticyclonic_Cell2(:,12); Cyclonic_Cell2(:,11); Cyclonic_Cell2(:,12)]));

            %Set original field
            origField = centerField;

            %Calculate Similarity Scores
            simTable = zeros(4,4);
            [simTable(:,:), simScore(:,:)] = aeceSimilarityMatrix(Label_anti2, Label_cyclo2,...
                Label_anti1, Label_cyclo1, Nanti2, Ncyclo2, Nanti, Ncyclo,...
                ratioTest, ratioKey, origField);
            FieldsSimilarity = [
                '(1,1) AE Number Test            '; ...
                '(1,2) CE Number Test            '; ...
                '(1,3) AE Number Key             '; ...
                '(1,4) CE Number Key             '; ...
                '(2,1) AE True Positive          '; ...
                '(2,2) AE False Positive         '; ...
                '(2,3) AE False Negative         '; ...
                '(2,4) AE True Negative          '; ...
                '(3,1) CE True Positive          '; ...
                '(3,2) CE False Positive         '; ...
                '(3,3) CE False Negative         '; ...
                '(3,4) CE True Negative          '; ...
                '(4,1) Number Error              '; ...
                '(4,2) Spatial Error Correct     '; ...
                '(4,3) Spatial Error Incorrect   '; ...
                '(4,4) Ratio Error               '];

            totalSims(fileNum,methCount) = simScore;
            totalRatio(fileNum,methCount) = simTable(4,4);
            totalSpatialC(fileNum,methCount) = simTable(4,2);
            totalSpatialR(fileNum,methCount) = simTable(4,3);
            totalNum(fileNum,methCount) = simTable(4,1);

            %Save Similarity Scores back into file
            save(filename_out,'simTable','simScore','FieldsSimilarity','-append');
            toc
        end %END LISTFOR
        methCount = methCount + 1;
    end %END OPTFOR

    %Display winnter messages
    [MaxVal,MaxInd] = max(mean(totalSims,1));
    [MaxValr,MaxIndr] = min(mean(totalRatio,1));
    [MaxValc,MaxIndc] = min(mean(totalSpatialR,1));
    [MaxVals,MaxInds] = min(mean(totalSpatialC,1));
    [MaxValn,MaxIndn] = min(mean(totalNum,1));
    if optimizeMethod
        ratioMessage = ["The ratio winner is " + masterMethods(1,MaxIndr) + " Centers and " + masterMethods(2,MaxIndr) + " Edges"];
        numMessage = ["The number winner is " + masterMethods(1,MaxIndn) + " Centers and " + masterMethods(2,MaxIndn) + " Edges"];
        spatialMessage = ["The spatial correct winner is " + masterMethods(1,MaxInds) + " Centers and " + masterMethods(2,MaxInds) + " Edges"];
        spatialcMessage = ["The spatial incorrect winner is " + masterMethods(1,MaxIndc) + " Centers and " + masterMethods(2,MaxIndc) + " Edges"];
        winnerMessage = ["The overall winner is " + masterMethods(1,MaxInd) + " Centers and " + masterMethods(2,MaxInd) + " Edges"];
        xToPlot = 1:length(masterMethods);
    elseif optimizeSTD
        ratioMessage = ["The ratio winner is " + num2str(stdFactor(MaxIndr))];
        numMessage = ["The number winner is " + num2str(stdFactor(MaxIndn))];
        spatialMessage = ["The spatial correct winner is " + num2str(stdFactor(MaxInds))];
        spatialcMessage = ["The spatial incorrect winner is " + num2str(stdFactor(MaxIndc))];
        winnerMessage = ["The best overall STD Factor is " + num2str(stdFactor(MaxInd))];
        xToPlot = stdFactor;
        winner = stdFactor(MaxInd);
    elseif optimizeK
        ratioMessage = ["The ratio winner is " + num2str(K(MaxIndr))];
        numMessage = ["The number winner is " + num2str(K(MaxIndn))];
        spatialMessage = ["The spatial correct winner is " + num2str(K(MaxInds))];
        spatialcMessage = ["The spatial incorrect winner is " + num2str(K(MaxIndc))];
        winnerMessage = ["The best overall K Factor is " + num2str(K(MaxInd))];
        xToPlot = K;
        winner = K(MaxInd);
    elseif optimizeSmooth
        ratioMessage = ["The ratio winner is " + num2str(smoothingFactor(MaxIndr))];
        numMessage = ["The number winner is " + num2str(smoothingFactor(MaxIndn))];
        spatialMessage = ["The spatial correct winner is " + num2str(smoothingFactor(MaxInds))];
        spatialcMessage = ["The spatial incorrect winner is " + num2str(smoothingFactor(MaxIndc))];
        winnerMessage = ["The best overall Smoothing Factor is " + num2str(smoothingFactor(MaxInd))];
        xToPlot = smoothingFactor;
        winner = smoothingFactor(MaxInd);
    elseif optimizeDist
        ratioMessage = ["The ratio winner is " + num2str(minDist(MaxIndr))];
        numMessage = ["The number winner is " + num2str(minDist(MaxIndn))];
        spatialMessage = ["The spatial correct winner is " + num2str(minDist(MaxInds))];
        spatialcMessage = ["The spatial incorrect winner is " + num2str(minDist(MaxIndc))];
        winnerMessage = ["The best overall Distance is " + num2str(minDist(MaxInd))];
        xToPlot = minDist;
    elseif optimizeRatio
        ratioMessage = ["The ratio winner is " + num2str(minRatio(MaxIndr))];
        numMessage = ["The number winner is " + num2str(minRatio(MaxIndn))];
        spatialMessage = ["The spatial correct winner is " + num2str(minRatio(MaxInds))];
        spatialcMessage = ["The spatial incorrect winner is " + num2str(minRatio(MaxIndc))];
        winnerMessage = ["The best overall Ratio is " + num2str(minRatio(MaxInd))];
        xToPlot = minRatio;
    elseif optimizeAandB
        ratioMessage = ["The ratio winner is " + masterAB(1,MaxIndr) + " A and " + masterAB(2,MaxIndr) + " B"];
        numMessage = ["The number winner is " + masterAB(1,MaxIndn) + " A and " + masterAB(2,MaxIndn) + " B"];
        spatialMessage = ["The spatial correct winner is " + masterAB(1,MaxInds) + " A and " + masterAB(2,MaxInds) + " B"];
        spatialcMessage = ["The spatial incorrect winner is " + masterAB(1,MaxIndc) + " A and " + masterAB(2,MaxIndc) + " B"];
        winnerMessage = ["The overall winner is " + masterAB(1,MaxInd) + " A and " + masterAB(2,MaxInd) + " B"];
        xToPlot = 1:length(masterAB);
    elseif optimizeSmoothAndRatio
        ratioMessage = ["The ratio winner is " + masterSR(1,MaxIndr) + " Smoothing and " + masterSR(2,MaxIndr) + " Ratio"];
        numMessage = ["The number winner is " + masterSR(1,MaxIndr) + " Smoothing and " + masterSR(2,MaxIndr) + " Ratio"];
        spatialMessage = ["The spatial correct winner is " + masterSR(1,MaxIndr) + " Smoothing and " + masterSR(2,MaxIndr) + " Ratio"];
        spatialcMessage = ["The spatial incorrect winner is " + masterSR(1,MaxIndr) + " Smoothing and " + masterSR(2,MaxIndr) + " Ratio"];
        winnerMessage = ["The overall winner is " + masterSR(1,MaxIndr) + " Smoothing and " + masterSR(2,MaxIndr) + " Ratio"];
        xToPlot = 1:length(masterSR);
        winner = masterSR(1,MaxIndr);
        winner2 = masterSR(2,MaxIndr);
    end
    disp(ratioMessage);
    disp(numMessage);
    disp(spatialMessage);
    disp(spatialcMessage);
    disp(winnerMessage);
    disp(["Max Similarity Score: " + num2str(MaxVal)]);
    datetimeNow = datetime;
    filename_out = [basepath + "/SIMILARITIES_FINAL/" + "IO_" + titlestring + "_Comparisons.mat"];
    save(filename_out, "totalSims","totalRatio","totalSpatialC","totalSpatialR","totalNum",...
        "masterMethods","masterDirList","inputDir","datetimeNow","xToPlot","winner","winner2");
end %END RPTLOOP
