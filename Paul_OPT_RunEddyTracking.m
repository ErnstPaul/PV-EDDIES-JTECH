%% Paul_OPT_RunEddyTracking.m (version 1.2)
%Author: Paul Ernst
%Date Created: 6/2/2022
%Date of Last Update: 1/12/2023
%Update History:
%PE 1/12/2023 - Cleared in final revision
%PE 9/23/2022 - Cleaned for final release
%PE 6/2/2022 - Created
%--------------------------------------
%Purpose: Runs the Eddy Tracking for the Optimization Project PVEDDIES
%Inputs: See Below
%Outputs: Tracked eddies in a domain
%--------------------------------------

%% STEP 0: DEFAULTS
%Ensure compatibility with LoadDefaultSettings and DefaultData:
% clear;
% input_dir = '/Volumes/NEMO1_local/NEMOv3.1/'; %input for NEMO
% Paul_OPT_CreateDefaultData(input_dir); %Run at least once to get the DefaultData into your workspace
% disp("Initial Calibration completed.")

%% INPUTS
yearsMAT = ["2016","2017","2019"]; %enter all years formatted like ["2019", "2020"]-- DOUBLE QUOTES ESSENTIAL
monthsMAT = [1 2 5 6 7 8 9 11 12]; %months need to be doubles, e.g. [9,1]
daysMAT = [1:31]; %days need to be doubles, e.g. [1 14 28]
latlonbounds = [30, -10, 80, 40]; % [N, S, E, W] lat long boundaries
yearmetadata = [20160101, 20191231]; %[YYYYMMDD, YYYYMMDD] start and end dates of the data, only metadata
rhoBounds = [1000 1025.5; 1026.95 1027.4]; %first col: rhomin; second col: rhomax; surf and depth
outputDir = 'NEM_FINAL'; %name of the directory we're going to be using
% Reference Location for density rescaling profile
refLon = 72; refLat = 0;
centerMethod = "PV_ISO";
edgeMethod = "PV_ISO";
%Parameters for methods set by optimization (reference your figures for this)
K = 0.9; smthCenter = 3; stdCenter = 1;
smthEdge = 6; stdEdge = 1;
minRad = 25; minAmp = 0; minlife = 0; minDist = 25; minRatio = 0.12; %MinDist should be 25 pixels for accurate localmax search
% Depth levels for eddies (iso is corresponding to rhoBounds above)
dlev = 1;
dleviso = [1,2];

%% STEP 1: PREPARING FILES/CREATING DATA
depthlevel = 1; %we use the surface for optimization
load('DefaultData_OPT.mat') %See Paul_OPT_CreateDefaultSettings.m
Paul_OPT_LoadDefaultSettings(); %See corresponding function
output_dir_interim = strcat(basepath, ['/EXTRACTION/INDIVIDUAL_FILES/' outputDir]);
input_dir = strcat(output_dir_interim, [num2str(depthlevel) + "/"]);
inputDir = char(input_dir);
[input_dir] = Paul_OPT_FilePreparation(latlonbounds, yearmetadata, yearsMAT, monthsMAT, daysMAT, depthlevel, rhoBounds, outputDir,...
    refLon,refLat);

%% STEP 2: RUNNING EDDY ID
for i = 1:length(dleviso)
Paul_OPT_MethodPicker(inputDir,centerMethod,edgeMethod,smthCenter,smthEdge,stdCenter,stdEdge,K,minRatio,minAmp,minRad,minDist,dlev,dleviso(i))

%% TRACKING STEPS
addpath([basepath + "/EDDYTRACKING/"])
output_dir_interim = strcat(basepath, ['/EXTRACTION/EDDY_PROPERTIES/' outputDir]);
input_dir = strcat(output_dir_interim, [num2str(dleviso(i)) + "/"]);
propDir = char(input_dir);
Paul_OPT_Step3_SS(inputDir,propDir,dleviso(i)); %''

%% Step 4. Calculate statistics for preparation of trajectories.
Paul_OPT_Step4_SS(inputDir,propDir,dleviso(i)); %''
output_dir_interim = strcat(basepath, ['/EXTRACTION/EDDY_TRAJECTORIES/' outputDir]);
input_dir = strcat(output_dir_interim, [num2str(dleviso(i)) + "/"]);
trajDir = char(input_dir);

%% Step 5. Calculate trajectories.
Paul_OPT_Step5_SS(inputDir,propDir,trajDir,dleviso(i)); %''

%% Step 6. Merge trajectories into master.
Paul_OPT_Step6_SS(inputDir,propDir,trajDir,minlife,minAmp,minRad,dleviso(i)); %''
end
