%% Paul_OPT_CreateDefaultData.m (version 1.3)
%Author: Paul Ernst
%Date Created: 2/9/2022
%Date of Last Update: 1/12/2023
%Update History:
%PE 1/12/2023 - Cleared in final revision
%PE 9/23/2022 - Cleaned for final release
%PE 6/2/2022 - Updated for OPT folder
%PE 2/9/2022 - Created
%--------------------------------------
%Purpose: Creates comprehensive matrix with all vital input data for this project
%Inputs: The input directory for CalibrateNEMODataStructure
%Outputs: Outputs a file with all variables that are important to the project
%--------------------------------------
function Paul_OPT_CreateDefaultData(input_dir)
%% Clear workspace for saving
clc; close all; % we'll be saving the whole workspace in a second here, no clutter

%% Years & Months as strings
monthsS = ["01"; "02"; "03"; "04"; "05"; "06"; "07"; "08"; "09";...
    "10"; "11"; "12"];
yearsS =  ["1993"; "1994"; "1995"; "1996"; "1997"; "1998"; "1999"; "2000";...
    "2001"; "2002"; "2003"; "2004"; "2005"; "2006"; "2007"; "2008"; "2009";...
    "2010"; "2011"; "2012"; "2013"; "2014"; "2015"; "2016"; "2017"; "2018"; "2019"];

%% Years & Months as characters
monthsC = ['01'; '02'; '03'; '04'; '05'; '06'; '07'; '08'; '09';...
    '10'; '11'; '12'];
yearsC =  ['1993'; '1994'; '1995'; '1996'; '1997'; '1998'; '1999'; '2000';...
    '2001'; '2002'; '2003'; '2004'; '2005'; '2006'; '2007'; '2008'; '2009';...
    '2010'; '2011'; '2012'; '2013'; '2014'; '2015'; '2016'; '2017'; '2018'; '2019'];

%% NEMO3.3.1 Years & Months as strings/chars
monthsSNEM = ["01"; "02"; "03"; "04"; "05"; "06"; "07"; "08"; "09";...
    "10"; "11"; "12"];
yearsSNEM =  ["2016"; "2017"; "2018"; "2019"; "2020"; "2021"];
monthsCNEM = ['01'; '02'; '03'; '04'; '05'; '06'; '07'; '08'; '09';...
    '10'; '11'; '12'];
yearsCNEM =  ['2016'; '2017'; '2018'; '2019'; '2020'; '2021'];

%% NEMO Input Structure
input_final_NEM = Paul_OPT_CalibrateNEMODataStructure(input_dir);

%% NEMO Lats, Lons, & Depths (read from file)
latsNEM = double(ncread(input_final_NEM{1,1,1}, 'latitude'));
lonsNEM = double(ncread(input_final_NEM{1,1,1}, 'longitude'));
depthsNEM = transpose(double(ncread(input_final_NEM{1,1,1}, 'depth')));

%% Paths
basepath = pwd;
funcpath = strcat(pwd, '/FUNCTIONS/');
matpath = strcat(pwd, '/SIMILARITIES_FINAL/');
input_dir_NEM = input_dir;
clear input_dir

%% Save workspace
save("DefaultData_OPT.mat")
end