%% Paul_OPT_TrajFigureDiags
%Author: Yves Morel
%Date Created: 02/2020
%Date of Last Update: 1/12/2023
%Update History:
%PE 1/12/2023 - Cleared in final revision
%PE 7/25/2022 - Revamped into new context
%--------------------------------------
%Purpose: Creates the figure figure outputs for the particle tracking suite (Figures 15-17)
%Inputs: The All-matrix for the particles from the suite
%Outputs: Figures as in Assene et al., 2020
%--------------------------------------
%% INPUTS
load('DefaultData_OPT.mat') %See Paul_OPT_CreateDefaultSettings.m
Paul_OPT_LoadDefaultSettings(); %See corresponding function
partPath = [basepath '/PARTICLE_TRACKING/'];
addpath(partPath)

%String to use for the run
titlestring = 'NEM_RSW_CHAGOS_BACKWARD_HIGH';

%Temporal variables
date2start  =  '20190529'; %YYYYMMDD day to start on in single quotes
days2integrate  =  100; %number of days to integrate forward or backward

%Spatial variables
latlonbounds  =  [-2, -13, 88, 65]; %the N S E W bounds of the overall area to be studied (domain)
depth2use = 34;
npart = 500;
fn = 18;
centerLon = 71;
delin = 1; %delineation for Ri number to show

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

%sorting by east and west - literally just select indexvortnew2-3-4 for different high/low/med pv
%figures (with adjustments to cmaps and stuff)
indexvortnew3east = indexvortnew2(ismember(indexvortnew2,find(lonalltraj(Ninit,:)>centerLon)));
east2plot = indexvortnew3east;
indexvortnew3west = indexvortnew2(ismember(indexvortnew2,find(lonalltraj(Ninit,:)<=centerLon)));
west2plot = indexvortnew3west;

%% Create figure of all particles + PV + depth evolution + their places on a T/S Diagram / Rho/PV Diagram
figure('units', 'normalized', 'outerposition', [0 0.0652777777777778 0.57890625 0.917361111111111], 'Color', 'w'); %fullscreen, invisible figure
t = tiledlayout(3,2,'TileSpacing','compact');
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
terrainMapX = Xgrid(isnan(nansToBlock));%CHAGOS ARCHIPELAGO
terrainMapY = Ygrid(isnan(nansToBlock));%CHAGOS ARCHIPELAGO
darkGrey1  = [0.4 0.4 0.4];%CHAGOS ARCHIPELAGO

% EAST TRAJECTORY PLOT
ax = nexttile(1);
hold on
%plot lines
for ipart = 1:length(east2plot)
    Xc = lonalltraj(:,east2plot(ipart)); Yc = latalltraj(:,east2plot(ipart)); Ac = PValltraj(:,east2plot(ipart));
    Xc(end+1) = NaN; Yc(end+1) = NaN; Ac(end+1) = NaN;
    patch(Xc,Yc,Ac,'EdgeColor','Interp','LineWidth',3) %colored lines
end
%terrain
scatter(terrainMapX,terrainMapY,48,darkGrey1,'filled','s')
%Mappingthings
text(71,-4.5,"(a)",'color','k','fontsize',32)
ylim([-10 -4])
xlim([70 82])
caxis([-1.75e-5 0])
grid on
box on
xlabel("Longitude (ºE)")
ylabel("Latitude (ºS)")
title("Trajectories of Eastern High PV Particles (s^-^1)")
set(gca,'fontsize',fn);
set(gca,'linewidth',4);
colormap(ax, jet)
colorbar();

% EAST DENSITY/PV PLOT
ax = nexttile(3);
hold on
for ipart = 1:length(east2plot)
    Xc = fliplr(timetraj); Yc = densityalltraj(:,east2plot(ipart)); Ac = PValltraj(:,east2plot(ipart));
    scatter(Xc,Yc,[],Ac,'filled')
end
xlim([-1 101])
xticks([11 42 72])
xticklabels([convertCharsToStrings(datestr(Tend+11)),...
    convertCharsToStrings(datestr(Tend+42)),...
    convertCharsToStrings(datestr(Tend+72))])
ylim([1026.8 1027.4])
caxis([-1.75e-5 0])
text(17,1027.35,"(b)",'color','k','fontsize',32)
colormap(ax, jet)
xlabel("Date")
ylabel("Potential Density (kg m^-^3)")
title("PV Evolution of Eastern High PV Particles (s^-^1)")
grid on
box on
set(gca,'fontsize',fn);
set(gca,'linewidth',4);
colorbar()

% EAST DENSITY/RI PLOT
ax = nexttile(5);
for ipart = 1:length(east2plot)
    Xc = fliplr(timetraj); Yc = densityalltraj(:,east2plot(ipart)); Ac = Rialltraj(:,east2plot(ipart));
    determiner = Ac<delin;
    Xc(~determiner) = []; Yc(~determiner) = []; Ac(~determiner) = [];
    hold on
    scatter(Xc,Yc,[],Ac./delin,'filled')
end
xlim([-1 101])
xticks([11 42 72])
xticklabels([convertCharsToStrings(datestr(Tend+11)),...
    convertCharsToStrings(datestr(Tend+42)),...
    convertCharsToStrings(datestr(Tend+72))])
ylim([1026.8 1027.4])
caxis([0 1])
colormap(ax, copper)
text(10,1027.3,"(c)",'color','k','fontsize',32)
xlabel("Date")
ylabel("Potential Density (kg m^-^3)")
title("Ri Evolution of Eastern High PV Particles")
grid on
box on
set(gca,'fontsize',fn);
set(gca,'linewidth',4);
colorbar();



% WEST TRAJECTORY PLOT
ax = nexttile(2);
hold on
%plot lines
for ipart = 1:length(west2plot)
    Xc = lonalltraj(:,west2plot(ipart)); Yc = latalltraj(:,west2plot(ipart)); Ac = PValltraj(:,west2plot(ipart));
    Xc(end+1) = NaN; Yc(end+1) = NaN; Ac(end+1) = NaN;
    patch(Xc,Yc,Ac,'EdgeColor','Interp','LineWidth',3) %colored lines
end
%terrain
scatter(terrainMapX,terrainMapY,60,darkGrey1,'filled','s')
%Mappingthings
text(65,-5,"(d)",'color','k','fontsize',32)
ylim([-10 -4])
xlim([63 76])
caxis([-1.75e-5 0])
grid on
box on
xlabel("Longitude (ºE)")
ylabel("Latitude (ºS)")
title("Trajectories of Western High PV Particles (s^-^1)")
set(gca,'fontsize',fn);
set(gca,'linewidth',4);
colormap(ax, jet)
colorbar();

% WEST DENSITY/PV PLOT
ax = nexttile(4);
hold on
for ipart = 1:length(west2plot)
    Xc = fliplr(timetraj); Yc = densityalltraj(:,west2plot(ipart)); Ac = PValltraj(:,west2plot(ipart));
    Ac2 = Rialltraj(:,west2plot(ipart));
    scatter(Xc,Yc,[],Ac,'filled')
end
xlim([-1 101])
xticks([11 42 72])
xticklabels([convertCharsToStrings(datestr(Tend+11)),...
    convertCharsToStrings(datestr(Tend+42)),...
    convertCharsToStrings(datestr(Tend+72))])
ylim([1026.85 1027.2])
caxis([-1.75e-5 0])
colormap(ax, jet)
text(30,1027.145,"(e)",'color','k','fontsize',32)
xlabel("Date")
ylabel("Potential Density (kg m^-^3)")
title("PV Evolution of Western High PV Particles")
grid on
box on
set(gca,'fontsize',fn);
set(gca,'linewidth',4);
colorbar()

% WEST DENSITY/RI PLOT
ax = nexttile(6);
for ipart = 1:length(west2plot)
    Xc = fliplr(timetraj); Yc = densityalltraj(:,west2plot(ipart)); Ac = Rialltraj(:,west2plot(ipart));
    determiner = (Ac<delin);
    Xc(~determiner) = []; Yc(~determiner) = []; Ac(~determiner) = [];
    hold on
    scatter(Xc,Yc,[],Ac./delin,'filled')
end
xlim([-1 101])
xticks([11 42 72])
xticklabels([convertCharsToStrings(datestr(Tend+11)),...
    convertCharsToStrings(datestr(Tend+42)),...
    convertCharsToStrings(datestr(Tend+72))])
ylim([1026.85 1027.2])
caxis([0 1])
colormap(ax, copper)
text(10,1027.1,"(f)",'color','k','fontsize',32)
xlabel("Date")
ylabel("Potential Density (kg m^-^3)")
title("Ri Evolution of Western High PV Particles")
grid on
box on
set(gca,'fontsize',fn);
set(gca,'linewidth',4);
colorbar();

%% save
filenamestring = ([basepath + "/FIGURES/FINAL_RENDER/" + titlestring + ".tiff"]);
filename2save = char(filenamestring);
%we use export_fig because it automatically crops and is nice
export_fig(filename2save,'-m1.5','-a4','-opengl','-nocrop'); %saves to local directory as PNG
