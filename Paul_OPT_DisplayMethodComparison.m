%% Paul_OPT_DisplayMethodComparison.m (version 1.2)
%Author: Paul Ernst
%Date Created: 6/10/2022
%Date of Last Update: 1/12/2023
%Update History:
%PE 1/12/2023 - Cleared in final revision
%PE 9/23/2022 - Cleaned up for final release
%PE 6/10/2022 - Created
%--------------------------------------
%Purpose: Displays 6 panels, each with a method from the optimized matrix (i.e. Figures 6-9)
%Inputs: Takes from COMPARISON_CENTERMETHOD_EDGEMETHOD folders after optimizemthod is run
%Outputs: 6 panels with centers and contours: methods chosen in the Plots variable below
%Note that this is for a specific plot so you must change things as you want different results
%% INPUTS
date2use = "20160101"; %YYYYMMDD
Plots = ["KEYFILES", "OW_T_OW_WA", "VORT_WA_VORT_WA", "LNAM_OW_WA", "VORT_WA_PV_ISO", "PV_ISO_PV_ISO"];
field2use = ["SSH", "OW", "VORT", "OW", "PV", "PV"];
latlonbounds = [13, 5, 80, 64]; %NSEW bounds
titlestring = ["OptimizedMethodMap" + date2use]; %title of figure
qs = 4; %quiverspacing for arrows
load("whiteorangecmap.mat")
annolat = 12; %lat/long for annotation
annolon = 76;
datfile = "/EXTRACTION/INDIVIDUAL_FILES/NEM_FULL1/pv_uv_" + date2use + ".mat";
fn = 28; %fontsize

%% PLOT PREPARATION
load('DefaultData_OPT.mat') %See Paul_OPT_CreateDefaultSettings.m
Paul_OPT_LoadDefaultSettings(); %See corresponding function

%Get our paths togethers
numpanels = length(Plots);
datename = ["comp_" + date2use + ".mat"];
pathname = strings(1,numpanels);
for i = 1:length(Plots)
    pathname(i) = [basepath + "/COMPARISON/" + Plots(i) + "/" + datename];
end
NL = latlonbounds(1);   SL = latlonbounds(2); WL = latlonbounds(4);   EL = latlonbounds(3);

%Grid
load([basepath + datfile], "X", "Y", "U", "V", "F_ISO");
Xgrid = zeros(length(X),length(Y));
Ygrid = zeros(length(X),length(Y));
for p = 1:length(Y)
    Xgrid(:,p) = X(:,1);
end
for p = 1:length(X)
    Ygrid(p,:) = Y(:,1);
end

%% Plot our results

panelstr = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)"];

%Get our tiledlayout up and running
figure('units', 'normalized', 'outerposition', [0 0.0652777777777778 0.6015625 0.917361111111111], 'Color', 'w');
t = tiledlayout(3,2,'TileSpacing','tight');

for nTile = 1:numpanels
    %next panel
    ax = nexttile();
    %load next panel
    load(pathname(nTile))

    %Reprocess centers
    count = 1;
    for i = 1:length(centerLabels{1,1})
        thisCenter = centerLabels{1,1}(count);
        [a, b] = ind2sub([length(X),length(Y)],thisCenter);
        thisX = X(a); thisY = Y(b);
        if ((thisY < SL) || (thisY > NL) || (thisX < WL) || (thisX > EL))
            centerLabels{1,1}(count) = [];
            count = count - 1;
        end
        count = count + 1;
    end
    count = 1;
    for i = 1:length(centerLabels{1,2})
        thisCenter = centerLabels{1,2}(count);
        [a, b] = ind2sub([length(X),length(Y)],thisCenter);
        thisX = X(a); thisY = Y(b);
        if ((thisY < SL) || (thisY > NL) || (thisX < WL) || (thisX > EL))
            centerLabels{1,2}(count) = [];
            count = count - 1;
        end
        count = count + 1;
    end
    %create map
    m_proj('mercator','longitude',[WL EL],'latitude',[SL NL]) %initialize map
    hold on
    %acquire plottable contours
    logiAnti = (Label_anti>0); logiCyclo = (Label_cyclo>0).*1;
    logiTotal = double(logiAnti)+logiCyclo;

    %HANDLE FIELDS - COLORMAPS
    switch field2use(nTile)
        case "PV"
            edgeField{1,1} = edgeField{1,1}-F_ISO; %PV PLOT
            rbflag = 1;
            owflag = 0;
        case "SSH"
            edgeField = cell(1,1);
            edgeField{1,1} = ZOS;
            rbflag = 0;
            owflag = 0;
        case "VORT"
            rbflag = 1;
            owflag = 0;
        case "OW"
            rbflag = 1;
            edgeField{1,1} = centerField{1,1} * -1; %negative vals
            owflag = 1;
    end

    %PCOLOR
    m_pcolor(Xgrid,Ygrid,edgeField{1,1}) %Look at the edge field as pcolor
    if rbflag && owflag
        m_scatter(Xgrid(centerLabels{1,1}),Ygrid(centerLabels{1,1}), 'xg', 'SizeData', 200, 'LineWidth', 3); %centers, AE
        m_scatter(Xgrid(centerLabels{1,2}),Ygrid(centerLabels{1,2}), '+g', 'SizeData', 200, 'LineWidth', 3); %centers, CE
    elseif rbflag
        m_scatter(Xgrid(centerLabels{1,1}),Ygrid(centerLabels{1,1}), 'xg', 'SizeData', 200, 'LineWidth', 3); %centers, AE
        m_scatter(Xgrid(centerLabels{1,2}),Ygrid(centerLabels{1,2}), '+g', 'SizeData', 200, 'LineWidth', 3); %centers, CE
    else
        m_scatter(Xgrid(centerLabels{1,1}),Ygrid(centerLabels{1,1}), '+k', 'SizeData', 200, 'LineWidth', 3); %centers, AE
        m_scatter(Xgrid(centerLabels{1,2}),Ygrid(centerLabels{1,2}), 'xk', 'SizeData', 200, 'LineWidth', 3); %centers, CE
    end

    %HANDLE CAXIS
    switch field2use(nTile)
        case "PV"
            caxis([ -2e-5 2e-5]);
            colormap(ax, redblue())
        case "OW"
            caxis([ -5e-10 0  ]);
            colormap(ax, flipud(CustomColormap))
        case "VORT"
            caxis([ -2e-5 2e-5]);
            colormap(ax, redblue())
        case "SSH"
            caxis([0 .75]) %Color limits for SSH
            colormap(ax, "jet")
    end
    shading interp
    m_coast('patch',[.5 .5 .5]);

    %HANDLE GRIDS
    switch nTile
        case 1
            title("SSH Method (m)")
                diffX = EL-WL; diffY = NL-SL;
                % Unit vector of 1 m/s
                [hpv5, htv5] = m_vec(100, EL-diffY*0.1, NL+diffX*0.1, 100, 0, 'k', 'key', '1 m s^{-1}',...
                    'shaftwidth', 4, 'headwidth', 20, 'headlength', 24);
                set(htv5,'FontSize',28);
        case 2
            title("OW_T Centers, OW_W_A Edges (s^-^1)")
        case 3
            title("VORT_W_A Centers, VORT_W_A Edges (s^-^1)")
        case 4
            title("LNAM Centers, OW_W_A Edges (s^-^1)")
        case 5
            title("VORT_W_A Centers, PV_I_S_O Edges (s^-^1)")
        case 6
            title("PV_I_S_O Centers, PV_I_S_O Edges (s^-^1)")
    end
    m_grid('box', 'fancy','fontsize',fn);

    %HANDLE ANNOTATION AND TEXT
    m_text(annolon,annolat,panelstr(nTile), 'Color', 'w','fontsize',36)
    set(gca,'fontsize',fn);
    set(gca,'LineWidth',2);

    %HANDLE CONTOURS
    m_contour(Xgrid,Ygrid,logiTotal, 'Color', 'k') %contours of eddies
    colorbar;

    %HANDLE QUIVERS
    h1 = m_quiver(Xgrid(1:qs:end,1:qs:end),Ygrid(1:qs:end,1:qs:end),U(1:qs:end,1:qs:end),V(1:qs:end,1:qs:end), 'k');
    set(h1,'AutoScale','on', 'AutoScaleFactor', 1.2)
    set(h1,'LineWidth',1.2)
end

%% Save and export
filenamestring = ([basepath + "/FIGURES/FINAL_RENDER/" + titlestring + ".tiff"]);
filename2save = char(filenamestring);
%we use export_fig because it automatically crops and is nice
export_fig(filename2save,'-m1.5','-a4','-opengl','-nocrop'); %saves to local directory as PNG

