%% Paul_OPT_3DColumn.m (version 2.1)
%Author: Paul Ernst
%Date Created: 6/21/2022
%Date of Last Update: 1/12/2023
%Update History:
%PE 1/12/2023 - Cleared in final revision
%PE 7/21/2022 - Finalized usage of geometry
%PE 6/21/2022 - Created
%--------------------------------------
%Purpose: Creates a 3d columnar plot of density, salinity, and PV with U/V at a certain point
%Inputs: The lat/lon and date to be grabbed for this purpose, as well as the file containing PV
%Outputs: 3d Column plot with 3 panels in the style of Sun et al., 2022, Fig. 10.
%--------------------------------------
%% Inputs
titlestring = "3DColumnChagos"; %title of resulting figure
yearsMAT = ["2019"]; %year to grab
monthsMAT = [5]; %month to grab
daysMAT = [29]; %day to grab
latlonbounds = [-5, -10, 76, 71.5]; %the N S E W bounds of the area in question
depthsList = [2 29 35 38 40 42]; %the depths to be surveyed here in model levels
qs = 3; %spacing of the quivers
runName = "NEM_FINAL1"; %name of the PV input file to be used here
zlims = [-2600 0]; %depth limits
caxes = {[1025.5 1027.8], [34.5 34.9], [-2e-5 2e-5]}; %colorbar axes
titles = {"Potential Density (kg m^-^3)", "Salinity (psu)", "Potential Vorticity Anomaly (s^-^1)"};
K = 0.9;
rhoBounds = [1026.95, 1027.4];
%Eddy characteristics, contours obtained from this file of the relevant date (careful on cyclo/anti
load('/Volumes/Lacie-SAN/SAN2/Paul_Eddies/Subsurface_Optimization/EXTRACTION/INDIVIDUAL_FILES/NEM_FINAL1/pv_uv_20190529.mat','Label_anti')
eddyNum = 3; %contour number assigned in the above file - obtain from cross referencing x/y
centroidX = 74; %center x/y
centroidY = -8; 

%% Load Data
Label_anti2use = Label_anti==eddyNum;
load('DefaultData_OPT.mat') %See Paul_OPT_CreateDefaultSettings.m
Paul_OPT_LoadDefaultSettings(); %See corresponding function
NL = latlonbounds(1);   SL = latlonbounds(2); WL = latlonbounds(4);   EL = latlonbounds(3);
NLnem = nearestpoint(NL, latsNEM);   SLnem = nearestpoint(SL, latsNEM);
ELnem = nearestpoint(EL, lonsNEM);   WLnem = nearestpoint(WL, lonsNEM);
depthlevel = 1; dZ = 50;
%Calculate RHO and grab salinity from corresponding file
yearMark = 1;
year2use = find(yearsMAT(yearMark)==yearsSNEM);
monthNum = 1;
month2use = monthsMAT(monthNum);
dayNum = 1;
day2useNEM = daysMAT(dayNum);
% What is our single file that we want to load right now?
file2use = input_final_NEM{day2useNEM,month2use,year2use};
% Load relevant variables from the file (T, S, U, V, Z)
% Reads only the specified section of file and the depth levels above/below for PV
% gradient calculations
ptemps = ncread(file2use, 'thetao',...
    [WLnem,SLnem,depthlevel,1], [ELnem-WLnem+1, NLnem-SLnem+1, dZ, 1]);
sals = ncread(file2use, 'so',...
    [WLnem,SLnem,depthlevel,1], [ELnem-WLnem+1, NLnem-SLnem+1, dZ, 1]);
pNEMgrid = zeros((NLnem-SLnem+1),length(depthsNEM));
for latLoop = 1:(NLnem-SLnem+1)
    pNEMgrid(latLoop,:) = gsw_p_from_z(depthsNEM.*-1,latsNEM(SLnem+latLoop,1));
end
SA = sals; % Absolute Salinity conversion semi-optional, differences are often insignificant
% Convert temperature to conservative temperature across the entire latitude band
sizemat = size(ptemps);
ctNEM = nan(sizemat);
for latLoop = 1:(NLnem-SLnem+1)
    ctNEM(:,latLoop,:) = gsw_CT_from_t(squeeze(sals(:,latLoop,:)),squeeze(ptemps(:,latLoop,:)),pNEMgrid(latLoop,:));
end
% Calculate in-situ density across the entire longitude band
rhoNEM = nan(sizemat);
for lonLoop = 1:(ELnem-WLnem+1)
    rhoNEM(lonLoop,:,:) = gsw_rho(squeeze(sals(lonLoop,:,:)),squeeze(ctNEM(lonLoop,:,:)),0);
end
RHO = rhoNEM;
%Grab PV from corresponding file
input_file = strcat(basepath, ["/EXTRACTION/INDIVIDUAL_FILES/" + runName]);
year = yearsSNEM(year2use); month = monthsS(month2use);
if day2useNEM < 10
    dayS = insertBefore(num2str(day2useNEM),1,"0");
else
    dayS = num2str(day2useNEM);
end
datecat = [year + month + dayS];
fileName = ["/pv_uv_" + datecat + ".mat"];
%tighten boundaries to those of the above
load([input_file + fileName], "PV", "U", "V", "X", "Y", "F_ISO");
NL2 = nearestpoint(NL, Y);   SL2 = nearestpoint(SL, Y);
EL2 = nearestpoint(EL, X);   WL2 = nearestpoint(WL, X);
F = F_ISO(WL2:EL2,SL2:NL2);
PV = PV(WL2:EL2,SL2:NL2,:);
U = U(WL2:EL2,SL2:NL2,:);
V = V(WL2:EL2,SL2:NL2,:);
X = X(WL2:EL2);
Y = Y(SL2:NL2);
Label_anti2use = Label_anti2use(WL2:EL2,SL2:NL2);

%% Grab isopycnal depths
rhomin = rhoBounds(1,1);
rhomax = rhoBounds(1,2);
%loop thru i, j (X, Y)
Zmax = nan(sizemat(1)-1,sizemat(2)-1);
Zmin = nan(sizemat(1)-1,sizemat(2)-1);
for i=1:sizemat(1)-1
    for j=1:sizemat(2)-1
        % Find z for rhomin
        if min(isnan(RHO(i,j,:)))==1
            Zmin(i,j)=NaN;
            Zmax(i,j)=NaN;
        else
            k0=min(find((RHO(i,j,:)-rhomin).^2==nanmin((RHO(i,j,:)-rhomin).^2)));
            % SECURITY TO DEAL WITH SURFACE
            % INVERSION OF DENSITY ...
            if rhomin<RHO(i,j,1)
                k0=1;
            end
            if RHO(i,j,k0)<rhomin
                kmin=k0;
                kp1=min(dZ,kmin+1);
            else
                if k0==1
                    kmin=1;
                    kp1=1;
                else
                    kmin=max(1,k0-1);
                    kp1=min(dZ,kmin+1);
                end
            end
            if kp1>kmin
                zmin=depthsNEM(kmin)+(rhomin-RHO(i,j,kmin))*(depthsNEM(kp1)-depthsNEM(kmin))/...
                    (RHO(i,j,kp1)-RHO(i,j,kmin));
            else
                zmin=depthsNEM(kmin);
            end
            %
            % Find z for rhomax
            k0=min(find((RHO(i,j,:)-rhomax).^2==nanmin((RHO(i,j,:)-rhomax).^2)));
            if RHO(i,j,k0)<rhomax
                if k0==dZ
                    kmax=dZ;
                    km1=kmax;
                else
                    kmax=min(dZ,k0+1);
                    km1=max(1,kmax-1);
                end
            else
                kmax=k0;
                km1=max(1,kmax-1);
            end
            if km1<kmax
                zmax=depthsNEM(kmax)+(rhomax-RHO(i,j,kmax))*(depthsNEM(km1)-depthsNEM(kmax))/...
                    (RHO(i,j,km1)-RHO(i,j,kmax));
            else
                zmax=depthsNEM(kmax);
            end
            % layer thickness
            Zmin(i,j)=zmin;
            Zmax(i,j)=zmax;
        end
    end
end

%% Construct grid
Xgrid = zeros(length(X),length(Y));
Ygrid = zeros(length(X),length(Y));
for p = 1:length(Y)
    Xgrid(:,p) = X(:,1);
end
for p = 1:length(X)
    Ygrid(p,:) = Y(:,1);
end

%% Create figure
figure('units', 'normalized', 'outerposition', [0 0 1 1], 'Color', 'w');
t = tiledlayout(1,3,'TileSpacing','tight');
datmat = {RHO, SA, PV};
%loop thru tiles
for nTile = 1:3
    %Create first depth layer
    ax = nexttile();
    dat2use = datmat{nTile}; %ensure data is correct
    hold on
    %loop thru all depths to be used here
    for nDepth = 1:length(depthsList)
        %Increment zdata
        thisDepth = depthsNEM(depthsList(nDepth));
        if nTile == 3
            colormap(ax, redblue) %REDBLUE
            dat2use(:,:,depthsList(nDepth)) = dat2use(:,:,depthsList(nDepth)) - F;
        else
            colormap(ax, jet)
        end
        h = pcolor(Xgrid,Ygrid,dat2use(:,:,depthsList(nDepth)));
        h.ZData = ones(size(h.ZData)) * (-1*thisDepth);
        modX = Xgrid(1:qs:end,1:qs:end);
        modY = Ygrid(1:qs:end,1:qs:end);
        modZ = ones(size(modX)) * (-1*thisDepth);
        modU = U(1:qs:end,1:qs:end,depthsList(nDepth));
        modV = V(1:qs:end,1:qs:end,depthsList(nDepth));
        modW = zeros(size(modU));
        q = quiver3(modX,modY,modZ,...
            modU,modV,modW,...
            'Color', 'k', 'linewidth', 1.1);
        shading interp
        caxis(caxes{nTile});
        colorbar
    end
    title(titles{nTile})
    zlim(zlims)
    %Ensure view and boxes are correctly calibrated
    view(gca,[-13 10]); %angle of viewed figure
    set(gca,'fontsize',22);
    set(gca,'LineWidth',4);
    box on
    grid on
    xlabel("Longitude (ºE)")
    ylabel("Latitude (ºN)") %change here if in S/N
    zlabel("Depth (m)")
    % Eddy plotting (mesh grid)
    if nTile == 3
        X2Edd = Xgrid(Label_anti2use>0);
        Y2Edd = Ygrid(Label_anti2use>0);
        middleDepth = nanmean(nanmean((Zmin+Zmax)./2*-1))*ones(size(Xgrid));
        middleDepth2Edd = middleDepth(Label_anti2use>0);
        meshLabel_cyclo2use = double(Label_anti2use);
        meshLabel_cyclo2use(meshLabel_cyclo2use==0) = NaN;
        mesh(Xgrid,Ygrid,meshLabel_cyclo2use.*middleDepth,'EdgeColor', 'k');
        [~,c] = contour(Xgrid,Ygrid,Label_anti2use,'EdgeColor', 'k');
        c.ContourZLevel = (nanmean(nanmean(middleDepth)));
        thisX = nearestpoint(centroidX,X);
        thisY = nearestpoint(centroidY,Y);
        depthsToCovermain = Zmin(thisX,thisY):1:Zmax(thisX,thisY);
        depthsToCoverup = 0:30:max(depthsToCovermain);
        depthsToCoverdown = min(depthsToCovermain):30:zlims(1)*-1;
        %solid line in the middle
        centerX2covermain = ones(1,length(depthsToCovermain))*centroidX;
        centerY2covermain = ones(1,length(depthsToCovermain))*centroidY;
        scatter3(centerX2covermain,centerY2covermain,-1*depthsToCovermain,100,'ok','filled')
        %dotted line above
        centerX2coverup = ones(1,length(depthsToCoverup))*centroidX;
        centerY2coverup = ones(1,length(depthsToCoverup))*centroidY;
        scatter3(centerX2coverup,centerY2coverup,-1*depthsToCoverup,100,'vk','filled')
        %dotted line below
        centerX2coverdown = ones(1,length(depthsToCoverdown))*centroidX;
        centerY2coverdown = ones(1,length(depthsToCoverdown))*centroidY;
        scatter3(centerX2coverdown,centerY2coverdown,-1*depthsToCoverdown,100,'^k','filled')
    end
    if nTile == 1
        modX2 = Xgrid(1:end-1,1:end-1); modY2 = Ygrid(1:end-1,1:end-1);
        mesh(modX2,modY2,-1*Zmin)
        mesh(modX2,modY2,-1*Zmax)
    end
end

%% Save Diagram
filenamestring = ([basepath + "/FIGURES/FINAL_RENDER/" + titlestring + ".tiff"]);
filename2save = char(filenamestring);
%we use export_fig because it automatically crops and is nice
export_fig(filename2save,'-nocrop');
