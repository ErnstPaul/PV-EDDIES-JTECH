%% Paul_OPT_LoadDefaultSettings.m (version 1.2)
%Author: Paul Ernst
%Date Created: 2/9/2022
%Date of Last Update: 1/12/2023
%Update History:
%PE 1/12/2023 - Cleared in final revision
%PE 6/2/2022 - Changed to OPT
%PE 2/9/2022 - Created
%--------------------------------------
%Purpose: Adds all the paths and configures all the settings that we need for this project
%Inputs: DefaultData file for this project
%Outputs: None
%--------------------------------------
function Paul_OPT_LoadDefaultSettings()
warning off;
% close all; clc; %you can uncomment these if you like
load('DefaultData_OPT.mat');
addpath(funcpath);
addpath(matpath);
addpath([basepath + "/DATA/"])
addpath(strcat(funcpath, 'm_map/private/')); 
addpath(strcat(funcpath, 'm_map/'))
addpath(strcat(funcpath, 'export_figs/')) 
addpath(strcat(funcpath, 'GSW/'))
addpath(strcat(funcpath, 'GSW/html/'))
addpath(strcat(funcpath, 'GSW/library/'))
addpath(strcat(funcpath, 'GSW/pdf')) 
addpath(strcat(funcpath, 'GSW/thermodynamics_from_t')) 
end

