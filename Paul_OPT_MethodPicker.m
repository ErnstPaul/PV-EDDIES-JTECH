%% Paul_OPT_MethodPicker.m (version 3.2)
%Author: Paul Ernst
%Date Created: 5/10/2022
%Date of Last Update: 1/12/2023
%Update History:
%PE 1/12/2023 - Cleared in final revision
%PE 9/23/2022 - Cleaned for final release
%PE 6/2/2022 - Adjusted to work for EdgeParams
%PE 5/24/2022 - Adjusted for PV_ISO corrections from Yves Morel
%PE 5/11/2022 - Added SSH detection and altered params
%PE 5/10/2022 - Created
%--------------------------------------
%Purpose: Takes a request for comparing eddy tracking methods and does so - OPERATIONALLY
%Inputs: All defaults, individual files created by file preparation steps
%        The contour/center methods we are trying are:
%           0: SSH (winding angle, ground truth method) (CALLED SSH)
%           1: Isopycnal PV (CALLED PV_ISO, smoothed by factor smoothingFactor)
%           2: Vorticity (CALLED VORT_T and VORT_WA, smoothed by factor smoothingFactor with stdFactor)
%           3: OW (raw) (CALLED OW_T and OW_WA, smoothed by factor smoothingFactor with stdFactor)
%           4: Local Normalized Angular Momentum (CALLED LNAM) (ONLY USABLE FOR CENTERS!)
%Outputs: Comparison table and separate files for each algorithm
%--------------------------------------
%% Inputs
function Paul_OPT_MethodPicker(inputDir,centerMethod,edgeMethod,smthCenter,smthEdge,stdCenter,stdEdge,K,minRatio,minAmp,minRad,minDist,dlev,dleviso)
load('DefaultData_OPT.mat') %See Paul_OPT_CreateDefaultSettings.m
Paul_OPT_LoadDefaultSettings(); %See corresponding function

%% Enable file reading
%Filter for usable files in the input directory
list = dir([ inputDir '/' 'pv*.mat']);
configFile = dir([ inputDir '/' 'Config*.mat']);
load([configFile.folder '/' configFile.name]);

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

disp("Beginning eddy extraction: " + centerMethod + " Centers with " + edgeMethod + " Edges.")

%Configure list
for fileNum = 1:length(list)
    tic
    disp(['Extracting eddies for:    ' list(fileNum).name ''])
    filename = [inputDir  list(fileNum).name];
    %Load and reaffirm presence of all necessary variables
    load(filename)
    EKE = double(EKE_ISO(:,:,dleviso));
    U = double(U(:,:,dlev));
    V = double(V(:,:,dlev));
    U_ISO = double(U_ISO(:,:,dleviso));
    V_ISO = double(V_ISO(:,:,dleviso));
    PV_ISO = double(PV_ISO(:,:,dleviso));
    VORT = double(VORT(:,:,dlev));
    LOW = double(LOW(:,:,dlev));
    ZOS = double(ZOS(:,:));
    LNAM = double(LNAM(:,:,dlev));
    LOW_ISO = double(LOW_ISO(:,:,dleviso));
    LNAM_ISO = double(LNAM_ISO(:,:,dleviso));

    %Center detection switch cases
    switch centerMethod
        case "SSH"
            centerField = ZOS;
        case "PV_ISO"
            centerField = cell(1,2);
            centerField{1,1} = PV_ISO;
            centerField{1,2} = smthCenter;
        case "VORT_WA"
            centerField = cell(1,2);
            centerField{1,1} = VORT;
            centerField{1,2} = smthCenter;
        case "VORT_T"
            centerField = cell(1,4);
            centerField{1,1} = VORT;
            centerField{1,2} = VORTSTATS;
            centerField{1,3} = stdCenter;
            centerField{1,4} = smthCenter;
        case "OW_WA"
            centerField = cell(1,2);
            centerField{1,1} = OW;
            centerField{1,2} = smthCenter;
        case "OW_T"
            centerField = cell(1,4);
            centerField{1,1} = OW;
            centerField{1,2} = OWSTATS;
            centerField{1,3} = stdCenter;
            centerField{1,4} = smthCenter;
        case "LNAM"
            centerField = cell(1,5);
            centerField{1,1} = LOW;
            centerField{1,2} = LNAM;
            centerField{1,3} = U;
            centerField{1,4} = V;
            centerField{1,5} = K;
        case "LNAM_ISO"
            centerField = cell(1,5);
            centerField{1,1} = LOW_ISO;
            centerField{1,2} = LNAM_ISO;
            centerField{1,3} = U_ISO;
            centerField{1,4} = V_ISO;
            centerField{1,5} = K;
    end

    %Edge detection switch cases
    switch edgeMethod
        case "SSH"
            edgeField = ZOS;
        case "PV_ISO"
            edgeField = cell(1,2);
            edgeField{1,1} = PV_ISO;
            edgeField{1,2} = smthEdge;
        case "VORT_WA"
            edgeField = cell(1,2);
            edgeField{1,1} = VORT;
            edgeField{1,2} = smthEdge;
        case "VORT_T"
            edgeField = cell(1,4);
            edgeField{1,1} = VORT;
            edgeField{1,2} = VORTSTATS;
            edgeField{1,3} = stdEdge; %thisSTD
            edgeField{1,4} = smthEdge;
        case "OW_WA"
            edgeField = cell(1,3);
            edgeField{1,1} = OW;
            edgeField{1,2} = smthEdge;
            edgeField{1,3} = VORT; %to distinguish AE (neg) from CE (pos)
        case "OW_T"
            edgeField = cell(1,5);
            edgeField{1,1} = OW;
            edgeField{1,2} = OWSTATS;
            edgeField{1,3} = stdEdge; %thisSTD
            edgeField{1,4} = smthEdge;
            edgeField{1,5} = VORT; %to distinguish AE (neg) from CE (pos)
    end

    %Run algorithm that we've constructed
    [Anticyclonic_Cell,Cyclonic_Cell,Nanti,Ncyclo,Label_anti,Label_cyclo,centerLabels] = ...
        Paul_OPT_EddyExtractionMaster(X,Y,centerField,edgeField,minDist,...
        minAmp,minRad,minRatio,EKE,centerMethod,edgeMethod);
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
    save(filename,'Anticyclonic_Cell','Cyclonic_Cell','Nanti','Ncyclo','centerLabels',...
        'Label_anti','Label_cyclo','FieldsCells', '-append');
    toc
end %END LISTFOR
end %END FUNCTION

