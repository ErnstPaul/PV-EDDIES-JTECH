%% Paul_OPT_Step1_PT.m (version 2.1)
%Author: Yves Morel
%Date Created: 02/2020
%Date of Last Update: 1/12/2023
%Update History:
%PE 1/12/2023 - Cleared in final revision
%PE 9/23/2022 - Cleaned for final release
%PE 6/23/2022 - Generalized and revamped portability
%--------------------------------------
%Purpose: Runs vorticity calculation for a given particle
%Inputs: Inputs from Paul_OPT_RunParticleTracking
%Outputs: Vorticity for a given particle
%--------------------------------------

function [vorttraj] = Paul_OPT_Step1_PT(lontraj,lattraj,zedtraj,date2start)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate vorticity at particle position
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('DefaultData_OPT.mat') %See Paul_OPT_CreateDefaultSettings.m
% Paul_OPT_LoadDefaultSettings(); %See corresponding function

%find indices within which particle is postioned before step
icel = nearestpoint(lontraj,lonsNEM);
jcel = nearestpoint(lattraj,latsNEM);
kcel = nearestpoint(zedtraj,depthsNEM);

%Handling extreme edge cases at grid edge (should not happen most of the time)
if icel == 1
    icel = icel+1;
end
if jcel == 1
    icel = icel+1;
end
if kcel == 1
    icel = icel+1;
end

%grab the file to use
day2use = str2num(date2start(7:8));
month2use = str2num(date2start(5:6));
year2use = str2num(date2start(1:4));
minYear = 2015; %NEMO model v3.3.1
year2use = year2use - minYear;
file2use = input_final_NEM{day2use,month2use,year2use}; %file determination

%Load U and V
u3D = ncread(file2use,'uo',[icel-1,jcel-1,kcel-1,1],[3,3,3,1]);
v3D = ncread(file2use,'vo',[icel-1,jcel-1,kcel-1,1],[3,3,3,1]);

%Calculate indices and diffs in xs and ys
i=2; j=2; k=2;
dxv = ones(3,3)*(111.32/12)*1000*cosd(lattraj); %DXDIFF
dyu = ones(3,3)*(110.573/12)*1000; %DYDIFF

%Calculate vorticity and return
zetaz1=(v3D(i+1,j,k)-v3D(i,j,k))/dxv(i+1,j);
zetaz1=zetaz1-(u3D(i,j+1,k)-u3D(i,j,k))/dyu(i,j+1);
zetaz2=(v3D(i+1,j,k+1)-v3D(i,j,k+1))/dxv(i+1,j);
zetaz2=zetaz2-(u3D(i,j+1,k+1)-u3D(i,j,k+1))/dyu(i,j+1);
vorttraj = 0.5*(zetaz1+zetaz2);

end
