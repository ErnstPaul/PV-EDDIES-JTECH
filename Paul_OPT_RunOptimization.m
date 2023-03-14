%% Paul_OPT_RunOptimization.m (version 1.2)
%Author: Paul Ernst
%Date Created: 6/2/2022
%Date of Last Update: 1/12/2023
%Update History:
%PE 1/12/2023 - Cleared in final revision
%PE 9/23/2022 - Cleaned for final release
%PE 6/2/2022 - Created
%--------------------------------------
%Purpose: Runs the Optimization Project
%Inputs: See Below
%Outputs: Optimized Methodology for Subsurface tracking in a domain
%--------------------------------------
% STEPS:
% 0: Create Defaults and ensure paths
% 1: Create Data
% 2: Create Key
% 3: Center Params
% 4: Center Smoothing
% 5: Edge Params
% 6: Edge Smoothing + Ratio
% 7: Method Optimization
%% STEP 0: DEFAULTS
%Ensure compatibility with LoadDefaultSettings and DefaultData:
% clear;
% input_dir = '/Volumes/NEMO1_local/NEMOv3.1/'; %input for NEMO
% Paul_OPT_CreateDefaultData(input_dir); %Run at least once to get the DefaultData into your workspace
% disp("Initial Calibration completed.")

%% INPUTS
load('DefaultData_OPT.mat') %See Paul_OPT_CreateDefaultSettings.m
Paul_OPT_LoadDefaultSettings(); %See corresponding function
yearsMAT = ["2016","2017","2019"]; %enter all years formatted like ["2019", "2020"]-- DOUBLE QUOTES ESSENTIAL
monthsMAT = [1 2 5 6 7 8 9 11 12]; %months need to be doubles, i.e. [9,1]
daysMAT = [1:31]; %days need to be doubles, i.e. [1 14 28]
latlonbounds = [30, -10, 80, 40]; % [N, S, E, W] lat long boundaries
yearmetadata = [20160101, 20191231]; %[YYYYMMDD, YYYYMMDD] start and end dates of the data, only metadata
depthlevel = 1; %we use the surface for optimization
rhoBounds = [1000 1025.5]; %first col: rhomin; second col: rhomax
outputDir = 'NEM_FINAL1'; %name of the directory we're going to be using, single quotes essential
% Reference Location for density rescaling profile
refLon = 72; refLat = 0;

%% STEP 1: PREPARING FILES/CREATING DATA
output_dir_interim = strcat(basepath, ['/EXTRACTION/INDIVIDUAL_FILES/' outputDir]);
input_dir = strcat(output_dir_interim, [num2str(depthlevel) + "/"]);
input_dir = char(input_dir);
disp("Step 1: Preparing Files (calculating PV/LNAM/etc.)")
[input_dir] = Paul_OPT_FilePreparation(latlonbounds, yearmetadata, yearsMAT, monthsMAT, daysMAT, depthlevel, rhoBounds, outputDir, refLon, refLat);

%% STEP 2: CREATE KEY FILES FOR TRACKING
disp("Step 2: Creating Key Files.")
Paul_OPT_Key(input_dir);

%% STEP 3: CENTER PARAM OPTIMIZATION
disp("Step 3: Optimizing Center Parameters.")
%STD - OW_T and VORT_T
%K - LNAM
Paul_OPT_CenterParams(input_dir);
%Display and save results
Paul_OPT_DisplayOptimizedCenterParams();

%% STEP 4: CENTER SMOOTHING OPTIMIZATION
%Smoothing for all center methods
disp("Step 4: Optimizing Center Smoothing.")
Paul_OPT_CenterSmoothing(input_dir);
%Display and save results
Paul_OPT_DisplayOptimizedCenterSmoothing();

%% STEP 5: EDGE PARAM OPTIMIZATION
disp("Step 5: Optimizing Edge Parameters.")
%STD - OW_T and VORT_T
Paul_OPT_EdgeParams(input_dir);
%Display and save results
Paul_OPT_DisplayOptimizedEdgeParams();

%% STEP 6: EDGE SMOOTHING/RATIO OPTIMIZATION
%Smoothing/Ratios for all edge methods
disp("Step 6: Optimizing Edge Smoothing and Ratios.")
Paul_OPT_EdgeSmoothingAndRatio(input_dir);
%Display and save results
Paul_OPT_DisplayOptimizedEdgeSmoothAndRatio();

%% STEP 7: METHOD OPTIMIZATION
%Smoothing for all center methods
disp("Step 7: Optimizing Methods.")
Paul_OPT_MethodOptimization(input_dir);
%Display and save results
Paul_OPT_DisplayOptimizedMethodsFull();

%% DONE
disp("Subsurface Eddy Tracking Optimization Done. Shutting Down.")