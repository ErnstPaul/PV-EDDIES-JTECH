%% Paul_OPT_Key.m (version 3.2)
%Author: Paul Ernst
%Date Created: 5/10/2022
%Date of Last Update: 1/12/2023
%Update History:
%PE 1/12/2023 - Cleared in final revision
%PE 9/23/2022 - Cleaned up code for final release
%PE 6/2/2022 - Adjusted to work for Key Files specifically
%PE 5/24/2022 - Adjusted for PV_ISO corrections from Yves Morel
%PE 5/11/2022 - Added SSH detection and altered params
%PE 5/10/2022 - Created
%--------------------------------------
%Purpose: Takes a request for comparing eddy tracking methods and does so (PRODUCES KEY SSH FILES)
%Inputs: All defaults, individual files created by file preparation steps
%        The contour/center methods we are trying are:
%           0: SSH (winding angle, ground truth method) (CALLED SSH)
%           1: Isopycnal PV (CALLED PV_ISO, smoothed by factor smoothingFactor)
%           2: Vorticity (CALLED VORT_T and VORT_WA, smoothed by factor smoothingFactor with stdFactor)
%           3: OW (raw) (CALLED OW_T and OW_WA, smoothed by factor smoothingFactor with stdFactor)
%           4: Local Normalized Angular Momentum (CALLED LNAM) (ONLY USABLE FOR CENTERS!)
%Outputs: Comparison table and separate files for each algorithm
%% Inputs
function Paul_OPT_Key(inputDir)
load('DefaultData_OPT.mat') %See Paul_OPT_CreateDefaultSettings.m
Paul_OPT_LoadDefaultSettings(); %See corresponding function

centerMethods = "SSH"; edgeMethods = "SSH";

createPairs = false;
%Factor for thresholds
stdFactor = 1; %0 - 2
smoothingFactor = 0; %0 - 15
K = 0.85; %0 - 1, factor for LNAM (0.85 max)
minAmp = 0;
minRatio = 0;
minRad = 25; %minimum radius (we're looking at mesoscale only, so 25+)
minDist = 25;

%Logicals for if we're optimizing something, pick only one of the following:
optimizeSTD = false;
optimizeSmooth = false;
optimizeK = false;
optimizeMethod = true;
optimizeDist = false;
optimizeRatio = false;
optimizeSmoothAndRatio = false;

%INPUT FILES
% - SEE Paul_FSS_Step1_ComparisonsMO_SS.m - First name of "pv"
%These files will have:
%Isopycnally Averaged PV_ISO at Surface
%Isopycnally Averaged VORT at Surface
%Surface: U, V, S, T, RHO, O-W Parameter, SSH, Vorticity, LNAM, EKE, X and Y

%OUTPUT FILES
%Location we're sending each of our file lists out to
outputDir = strcat(basepath, '/COMPARISON/');
[status,message,messageid] = mkdir(outputDir);


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

for optN = 1:optLen
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
        EKE = double(EKE(:,:));
        U = double(U(:,:));
        V = double(V(:,:));
        PV_ISO = double(PV_ISO(:,:));
        VORT = double(VORT(:,:));
        OW = double(OW(:,:));
        LOW = double(LOW(:,:));
        ZOS = double(ZOS(:,:));

        %Center detection switch cases
                centerField = ZOS;

        %Edge detection switch cases
                edgeField = ZOS;
           
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
    end %END LISTFOR
    methCount = methCount + 1;
end %END OPTFOR
end