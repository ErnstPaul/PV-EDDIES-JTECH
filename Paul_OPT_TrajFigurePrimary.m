%% Paul_OPT_TrajFigurePrimary
%Author: Yves Morel
%Date Created: 02/2020
%Date of Last Update: 1/12/2023
%Update History:
%PE 1/12/2023 - Cleared in final revision
%PE 7/25/2022 - Revamped into new context
%--------------------------------------
%Purpose: Creates the primary figure outputs for the particle tracking suite (Figure 14)
%Inputs: The All-matrix for the particles from the suite
%Outputs: Figures as in Assene et al., 2020
%--------------------------------------
%% INPUTS
load('DefaultData_OPT.mat') %See Paul_OPT_CreateDefaultSettings.m
Paul_OPT_LoadDefaultSettings(); %See corresponding function
partPath = [basepath '/PARTICLE_TRACKING/'];
addpath(partPath)

%String to use for the run
titlestring = 'NEM_RSW_CHAGOS_BACKWARD2';

%Temporal variables
date2start  =  '20190529'; %YYYYMMDD day to start on in single quotes
days2integrate  =  100; %number of days to integrate forward or backward

%Spatial variables
latlonbounds  =  [-2, -13, 88, 65]; %the N S E W bounds of the overall area to be studied (domain)
depth2use = 34; %depth to use for things at eddy center depth
npart = 500; %number of particles
fn = 22;

%% Load the all_file and grab particle pvs
output_dir = [partPath 'RESULTS/' titlestring];
file4load = [partPath 'RESULTS/' titlestring '/ALL'];
load(file4load);
PValltraj = (PValltraj);
Ninit = length(timetraj);
maxpv0=1;
indexvort=find(max(abs(PValltraj),[],1)<= maxpv0 );
medPVvort=median(PValltraj(1,indexvort));
minPVvort=min(PValltraj(1,indexvort));
maxPVvort=max(PValltraj(1,indexvort));
dPV=maxPVvort-minPVvort;
deltapv=2.*sqrt(median(abs(PValltraj(1,indexvort)-medPVvort).^2));
for index=1:4
    if index==2
        symbol='+';
        indexvortnew2=find(PValltraj(Ninit,:)-medPVvort>= deltapv & max(abs(PValltraj),[],1)<= maxpv0);
        addtext=' particles with high initial PV';
        % Color/symbol = red +
    elseif index==3
        symbol='<';
        indexvortnew3=find(PValltraj(Ninit,:)-medPVvort<=-deltapv & max(abs(PValltraj),[],1)<= maxpv0);
        addtext=' particles with low initial PV';
        % Color/symbol = blue <
    else
        symbol='d';
        indexvortnew4=find(abs(PValltraj(Ninit,:)-medPVvort)<=deltapv & max(abs(PValltraj),[],1)<= maxpv0);
        addtext=' particles with small PV variations';
        % Color/symbol = green d
    end
end

%% Create figure of all particles + PV + depth evolution + their places on a T/S Diagram / Rho/PV Diagram
figure('units', 'normalized', 'outerposition', [0 0.0652777777777778 0.358203125 0.917361111111111], 'Color', 'w'); %fullscreen, invisible figure
t = tiledlayout(4,2,'TileSpacing','compact');

%Coordinate grid from bounds
NL = latlonbounds(1);   SL = latlonbounds(2); WL = latlonbounds(4);   EL = latlonbounds(3);
NLnem = nearestpoint(NL, latsNEM);   SLnem = nearestpoint(SL, latsNEM);
ELnem = nearestpoint(EL, lonsNEM);   WLnem = nearestpoint(WL, lonsNEM);
Y = latsNEM(SLnem:NLnem-1,1);
X = lonsNEM(WLnem:ELnem-1,1);
Xgrid = zeros(length(X),length(Y));
Ygrid = zeros(length(X),length(Y));
for p = 1:length(Y)
    Xgrid(:,p) = X(:,1);
end
for p = 1:length(X)
    Ygrid(p,:) = Y(:,1);
end
nansToBlock = ncread(input_final_NEM{1,1,1},'so',[WLnem,SLnem,depth2use,1], [ELnem-WLnem, NLnem-SLnem, 1, 1]);
terrainMapX = Xgrid(isnan(nansToBlock)); %CHAGOS ARCHIPELAGO
terrainMapY = Ygrid(isnan(nansToBlock)); %CHAGOS ARCHIPELAGO
darkGrey1  = [0.4 0.4 0.4]; %CHAGOS ARCHIPELAGO

% PARTICLE TRAJECTORY PANEL
nexttile([1 2])
hold on
%highPV
scatter(lonalltraj(Ninit,indexvortnew2),latalltraj(Ninit,indexvortnew2),'r+')
%midPV
scatter(lonalltraj(Ninit,indexvortnew4),latalltraj(Ninit,indexvortnew4),'gd')
%lowPV
scatter(lonalltraj(Ninit,indexvortnew3),latalltraj(Ninit,indexvortnew3),'b<')
%endingloc
scatter(lonalltraj(1,:),latalltraj(1,:),'k*')
%terrain
scatter(terrainMapX,terrainMapY,48,darkGrey1,'filled','s') %CHAGOS ARCHIPELAGO
%Mappingthings
text(66,-3.5,"(a)",'color','k','fontsize',32)
ylim([-10 -2])
xlim([64 82])
grid on
box on
xlabel("Longitude (ºE)")
ylabel("Latitude (ºS)")
title("Initial (color) and Final (black) Particle Positions")
set(gca,'fontsize',fn);
set(gca,'linewidth',4);

% T/S DIAGRAM
nexttile;
%initial 
scatter(salinityalltraj(Ninit,indexvortnew2),temperaturealltraj(Ninit,indexvortnew2),'r+')
hold on
scatter(salinityalltraj(Ninit,indexvortnew3),temperaturealltraj(Ninit,indexvortnew3),'b<')
scatter(salinityalltraj(Ninit,indexvortnew4),temperaturealltraj(Ninit,indexvortnew4),'gd')
%final
scatter(salinityalltraj(1,indexvortnew2),temperaturealltraj(1,indexvortnew2),'k+')
scatter(salinityalltraj(1,indexvortnew3),temperaturealltraj(1,indexvortnew3),'k<')
scatter(salinityalltraj(1,indexvortnew4),temperaturealltraj(1,indexvortnew4),'kd')
xlim([34.7 34.95])
ylim([5 10.5])
text(34.72,9.5,"(b)",'color','k','fontsize',32)
ylabel("Conservative Temperature (ºC)")
xlabel("Salinity (psu)")
grid on
box on
set(gca,'fontsize',fn-4);
set(gca,'linewidth',4);

% RHO/PV DIAGRAM
nexttile;
%initial 
scatter(PValltraj(Ninit,indexvortnew2),densityalltraj(Ninit,indexvortnew2),'r+')
hold on
scatter(PValltraj(Ninit,indexvortnew3),densityalltraj(Ninit,indexvortnew3),'b<')
scatter(PValltraj(Ninit,indexvortnew4),densityalltraj(Ninit,indexvortnew4),'gd')
%final
scatter(PValltraj(1,indexvortnew2),densityalltraj(1,indexvortnew2),'k+')
scatter(PValltraj(1,indexvortnew3),densityalltraj(1,indexvortnew3),'k<')
scatter(PValltraj(1,indexvortnew4),densityalltraj(1,indexvortnew4),'kd')
xlim([-2.75e-5 -0.25e-5])
ylim([1026.65, 1027.5])
text(-2.5e-5,1027.4,"(c)",'color','k','fontsize',32)
xlabel("PV (s^-^1)")
ylabel("Potential Density (kg m^-^3)")
grid on
box on
set(gca,'fontsize',fn-4);
set(gca,'linewidth',4);

% DEPTH EVOLUTION PANEL
nexttile([1 2])
for i = 1:length(indexvortnew2)
    scatter(fliplr(timetraj),-1*zedalltraj(:,indexvortnew2(i)),0.25,'r+')
    hold on
end
for i = 1:length(indexvortnew3)
    scatter(fliplr(timetraj),-1*zedalltraj(:,indexvortnew3(i)),0.25,'b<')
end
for i = 1:length(indexvortnew4)
    scatter(fliplr(timetraj),-1*zedalltraj(:,indexvortnew4(i)),0.25,'gd')
end
xlim([-1 101])
xticks([11 42 72])
xticklabels([convertCharsToStrings(datestr(Tend+11)),...
    convertCharsToStrings(datestr(Tend+42)),...
    convertCharsToStrings(datestr(Tend+72))])
ylim([-1100 -400])
text(93,-480,"(d)",'color','k','fontsize',32)
xlabel("Date")
ylabel("Depth (m)")
grid on
box on
set(gca,'fontsize',fn);
set(gca,'linewidth',4);

% PV EVOLUTION PANEL
nexttile([1 2])
for i = 1:length(indexvortnew2)
    scatter(fliplr(timetraj),PValltraj(:,indexvortnew2(i)),0.25,'r+')
    hold on
end
for i = 1:length(indexvortnew3)
    scatter(fliplr(timetraj),PValltraj(:,indexvortnew3(i)),0.25,'b<')
end
for i = 1:length(indexvortnew4)
    scatter(fliplr(timetraj),PValltraj(:,indexvortnew4(i)),0.25,'gd')
end
xlim([-1 101])
ylim([-3.5e-5 3e-5])
xticks([11 42 72])
text(90,2e-5,"(e)",'color','k','fontsize',32)
xticklabels([convertCharsToStrings(datestr(Tend+11)),...
    convertCharsToStrings(datestr(Tend+42)),...
    convertCharsToStrings(datestr(Tend+72))])
xlabel("Date")
ylabel("PV (s^-^1)")
grid on
box on
set(gca,'fontsize',fn);
set(gca,'linewidth',4);

% Save and export
filenamestring = ([basepath + "/FIGURES/FINAL_RENDER/" + titlestring + "1.tiff"]);
filename2save = char(filenamestring);
%we use export_fig because it automatically crops and is nice
export_fig(filename2save,'-m1.5','-a4','-opengl','-nocrop'); %saves to local directory as PNG
    






