%% Paul_OPT_TS_RSW.m (version 2.2)
%Author: Paul Ernst
%Date Created: 6/23/2022
%Date of Last Update: 1/12/2023
%Update History:
%PE 1/12/2023 - Cleared in final revision
%PE 9/23/2022 - Cleaned for final release
%PE 6/23/2022 - Fixed detrend fails
%PE 6/15/2022 - Created
%--------------------------------------
%Purpose: Creates a T/S diagram of the red sea water water mass + diagram of avg. top/bot density
%Inputs: latlonbounds, years, months, and days to use
%Outputs: T/S diagram of the water mass in question
%--------------------------------------
%% Inputs
load('DefaultData_OPT.mat') %See Paul_OPT_CreateDefaultSettings.m
Paul_OPT_LoadDefaultSettings(); %See corresponding function

titlestring = "RSW_GOA_TS_Diagram";

refday = 7; refmonth = 7; refyear = 3; refLon = 55; refLat = 13;

%Input bounds of time and space
yearsMAT = ["2016","2017","2019"]; %enter all years formatted like ["2019", "2020"]-- DOUBLE QUOTES ESSENTIAL
monthsMAT = [1 2 5 6 7 8 9 11 12]; %months need to be doubles [9,1]
daysMAT = [1:31]; %days need to be doubles [1 14 28]
latlonbounds = [15.2197, 10.3781, 51.2867, 42.5135]; % ditto but for calculated area
latlonbounds2 = [30, -10, 80, 40]; % ditto but for calculated area
lowbound = 1026; highbound = 1028;

%% Get data if needed
%T-S Diagram
[CT2Plot,SA2Plot] = Paul_OPT_TSDiagram(latlonbounds, yearsMAT, monthsMAT, daysMAT);
%save back just in case
save([titlestring + "TSPARAMS.mat"], "CT2Plot", "SA2Plot", "yearsMAT", "daysMAT", "monthsMAT", "latlonbounds");

%Algorithm
[maxRSW,topminRSW,botminRSW] = Paul_OPT_TSDiagramAlgor(latlonbounds2, yearsMAT, monthsMAT, daysMAT, lowbound, highbound);
%save back just in case
save([titlestring + "TSPARAMSALGOR.mat"], "maxRSW", "topminRSW", "botminRSW", "yearsMAT", "daysMAT", "monthsMAT", "latlonbounds2");

%% Loading Statements
load([titlestring + "TSPARAMS.mat"])
load([titlestring + "TSPARAMSALGOR.mat"])

%ridding garbage data and data too deep to count (check your depths)
CT2Plot(CT2Plot<=0) = NaN;
CT2Plot(:,44:end) = NaN;
SA2Plot(SA2Plot<=30) = NaN;
SA2Plot(:,44:end) = NaN;
%% Plot diagram

%Get diff S and diff T
ss = SA2Plot(:); tt = CT2Plot(:);
smin=min(ss)-0.01.*min(ss);
smax=max(ss)+0.01.*max(ss);
thetamin=min(tt)-0.1*max(tt);
thetamax=max(tt)+0.1*max(tt);
xdim=round((smax-smin)./0.1+1);
ydim=round((thetamax-thetamin)+1);
thetai=((1:ydim)-1)*1+thetamin;
si=((1:xdim)-1)*0.1+smin;

%calculate diff dens
dens=zeros(ydim,xdim);
for j=1:ydim
    for i=1:xdim
        dens(j,i) = gsw_rho(si(i),thetai(j),0);
    end
end

%caclulate diff spiciness
spic=zeros(ydim,xdim);
for j=1:ydim
    for i=1:xdim
        spic(j,i) = gsw_spiciness0(si(i),thetai(j));
    end
end

%% Create figure
figure('units', 'normalized', 'outerposition', [0 0 1 1], 'Color', 'w');
t = tiledlayout(1,2,'TileSpacing','tight');

% T/S Diagram contours
ax = nexttile();

xlabel('Salinity (psu)','FontWeight','bold','FontSize',32)
ylabel('Conservative Temperature (ºC)','FontWeight','bold','FontSize',32)

% plotting scatter plot of CT2Plot and SA2Plot
hold on;
s = scatter(ss,tt,'.r');
xlim([34 42.5]) %check this
ylim([5 34.5])
cont = [1000:0.5:1050];
cont([53 57]) = [];
[c,h]=contour(si,thetai,dens,cont,'--k','ShowText','on');
clabel(c,h,'LabelSpacing',1000,'FontSize',18);
[cc,hh] = contour(si,thetai,dens,[lowbound highbound], '-k', 'linewidth', 5,'ShowText','on');
clabel(cc,hh,'LabelSpacing',1000,'FontSize',24);
text(35.0995594713656,33.1563393708294, "(a)", 'FontSize',48) %check this

% median profile
medprofS = nanmedian(SA2Plot,1);
medprofT = nanmedian(CT2Plot,1);
plot(medprofS, medprofT, '-k', 'Linewidth', 9);
hold on

% standard deviation profiles +/-
stddevS = nanstd(SA2Plot,1);
stddevT = nanstd(CT2Plot,1);
stddevS(1,40:34) = 0.03;
stddevT(1,39:43) = 1;
plot(stddevS+medprofS,stddevT+medprofT, '-k', 'Linewidth', 3)
plot(medprofS-stddevS,medprofT-stddevT, '-k', 'Linewidth', 3)
set(gca,'linewi',3)
set(gca,'FontSize',32)
box on

%% Plot of min/max and spiciness ---------------------------
ax = nexttile;
xlabel('Salinity (psu)','FontWeight','bold','FontSize',32)
ylabel('Conservative Temperature (ºC)','FontWeight','bold','FontSize',32)

%DENSITY CONTOURS IN BLACK
hold on;
cont = [1000:0.5:1050];
[c,h]=contour(si,thetai,dens,cont,'--k','ShowText','on');
clabel(c,h,'LabelSpacing',1000,'FontSize',18);
xlim([34.5 36.5])
ylim([0.5 28])
%SPICINESS CONTOURS IN RED
[c,h]=contour(si,thetai,spic,[-50:0.5:50],'--r');
text(34.6306167400881, 25.8527168732126, "(b)", 'FontSize',48)

% PROFILE
refLonNEM = nearestpoint(refLon, lonsNEM); refLatNEM = nearestpoint(refLat, latsNEM);
% Grab reference profile data
file2use = input_final_NEM{refday,refmonth,refyear};
ptemps = ncread(file2use, 'thetao',...
    [refLonNEM,refLatNEM,1,1], [1, 1, 50, 1]);
SA2Plot2 = ncread(file2use, 'so',...
    [refLonNEM,refLatNEM,1,1], [1, 1, 50, 1]);
zRef = double(ncread(file2use, 'depth', 1, 50));
% Pressure from Depth
pNEMgrid = gsw_p_from_z(zRef.*-1,refLat);
% Conservative T
CT2Plot2 = gsw_CT_from_t(squeeze(SA2Plot2),squeeze(ptemps),pNEMgrid);
CT2Plot2 = squeeze(CT2Plot2); SA2Plot2 = squeeze(SA2Plot2);

%Plot profiles
plot(SA2Plot2, CT2Plot2, '-k', 'Linewidth', 3);
dens2plot = gsw_rho(SA2Plot2, CT2Plot2, 0);
densYes = (dens2plot >= lowbound) & (dens2plot <= highbound);
densYes(1:32) = 0; %depth cutoff
plot(SA2Plot2(densYes), CT2Plot2(densYes), '-k', 'Linewidth', 9);

% Trend Line
x = SA2Plot2; x(isnan(x)) = [];
y = CT2Plot2; y(isnan(y)) = [];
coefficients = polyfit(x, y, 1);
xFit = linspace(min(x), max(x), length(x));
yFit = polyval(coefficients , xFit);
hold on;
plot(xFit, yFit, ':b', 'LineWidth', 3);

% PLOT HIGH AND LOW POINTS
%Find spiciness profile
spic2plot = gsw_spiciness0(SA2Plot2,CT2Plot2);
%detrend
spic2find = detrend(spic2plot,'omitnan');
spic2find(~densYes) = NaN;
%find max
[maxSpic,maxSpicInd] = nanmax(spic2find);
%find min above max
low2find = spic2find;
low2find(maxSpicInd:end) = NaN;
[minSpic,minSpicInd] = nanmin(low2find);
%find min below max
low2find = spic2find;
low2find(1:maxSpicInd-1) = NaN;
[minSpic2,minSpicInd2] = nanmin(low2find);
%plot
scatter(SA2Plot2(maxSpicInd),CT2Plot2(maxSpicInd), 800, 'or', 'filled', 'LineWidth', 5)
scatter(SA2Plot2(minSpicInd),CT2Plot2(minSpicInd), 800, '+r', 'LineWidth', 10)
scatter(SA2Plot2(minSpicInd2),CT2Plot2(minSpicInd2), 800, '+r', 'LineWidth', 10)
set(gca,'linewi',3)
set(gca,'FontSize',32)
box on

disp(["Density bounds: " + num2str(dens2plot(minSpicInd)) + " " + num2str(dens2plot(maxSpicInd))])

%% Save Diagram
filenamestring = ([basepath + "/FIGURES/FINAL_RENDER/" + titlestring + ".tiff"]);
filename2save = char(filenamestring);
%we use export_fig because it automatically crops and is nice
export_fig(filename2save,'-m1.5','-a4','-opengl','-nocrop'); %saves to local directory as PNG

close all

%% Bottom left panel for mean maxes of water mass isopycnals
figure('units', 'normalized', 'outerposition', [0 0 1 1], 'Color', 'w');
t = tiledlayout(1,2,'TileSpacing','tight');
ax = nexttile;

%construct grid
load([titlestring + "TSPARAMSALGOR.mat"])
Xgrid = zeros(length(X),length(Y));
Ygrid = zeros(length(X),length(Y));
for p = 1:length(Y)
    Xgrid(:,p) = X(:,1);
end
for p = 1:length(X)
    Ygrid(p,:) = Y(:,1);
end

%plot map and pcolor
NL = latlonbounds2(1);   SL = latlonbounds2(2); WL = latlonbounds2(4);   EL = latlonbounds2(3);
m_proj('mercator','longitude',[WL EL],'latitude',[SL NL]) %initialize map
m_pcolor(Xgrid,Ygrid,maxRSW) %ADT color plot in cm
shading interp
hold on
m_coast('patch',[.5 .5 .5]);
m_grid('box', 'fancy','fontsize',40, 'xtick', [44 52 60 68 76], 'xticklabels', [44 52 60 68 76]);
m_text(41.5,28.5,"(a)",'color','w','fontsize',48);
title("RSW Maximum Densities (kg m^-^3)",'FontSize',40)
colormap jet;
colorbar('Fontsize', 30, 'Position', [0.469937018657365 0.0813204508856683 0.00896923134263494 0.842995169082126])
caxis([1026.95 1027.15])

% Bottom right panel for mean ranges between top and bottom water mass isopycnals

ax = nexttile;

%plot map and pcolor
m_proj('mercator','longitude',[WL EL],'latitude',[SL NL]) %initialize map
m_pcolor(Xgrid,Ygrid,abs(topminRSW-botminRSW)) %ADT color plot in cm
shading interp
hold on
m_coast('patch',[.5 .5 .5]);
m_grid('box', 'fancy','fontsize',40, 'yticklabels', [], 'xtick', [44 52 60 68 76], 'xticklabels', [44 52 60 68 76]);
m_text(41.5,28.5,"(b)",'color','w','fontsize',48);
title("RSW Density Range Width (kg m^-^3)",'FontSize',40)
colormap jet;
caxis([0.4 0.75])
colorbar('Fontsize', 30, 'Position', [0.930928394668154 0.0805152979066023 0.00852473033184631 0.845410628019324])

%% Save Diagram
filenamestring = ([basepath + "/FIGURES/FINAL_RENDER/" + titlestring + "2.tiff"]);
filename2save = char(filenamestring);
%we use export_fig because it automatically crops and is nice
export_fig(filename2save,'-m1.5','-a4','-opengl','-nocrop'); %saves to local directory as PNG