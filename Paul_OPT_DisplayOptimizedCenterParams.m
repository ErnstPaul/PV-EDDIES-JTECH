%% Paul_OPT_DisplayOptimizedCenterParams.m (version 1.1)
%Author: Paul Ernst
%Date Created: 5/25/2022
%Date of Last Update: 1/12/2023
%Update History:
%PE 1/12/2023 - Cleared in final revision
%PE 9/23/2022 - Cleaned for final release
%PE 5/25/2022 - Created
%--------------------------------------
%Purpose: Displays the optimized parameters of center methods without smoothing (Figure 1)
%Inputs: Plot sources from Comparison_Main
%Outputs: 3 panels with the following lines:
% -Center LNAM: K    -Center VORT_T: STD   -Center OW_T: STD, SMTH 
%% Inputs
titlestring = "OptimizedCenterParams";
PlotSources = ["OPT_C_VORT_T_STD", "OPT_C_OW_T_STD", "OPT_LNAM_K"];
xLabels = ["STD Factor", "STD Factor", "K"];

%% Plots
load('DefaultData_OPT.mat') %See Paul_OPT_CreateDefaultSettings.m
Paul_OPT_LoadDefaultSettings(); %See corresponding function

numPanels = length(PlotSources);
figure('units', 'normalized', 'outerposition', [0 0 1 1], 'Color', 'w'); %fullscreen, invisible figure
t = tiledlayout(1,numPanels,'TileSpacing','tight');

panelLabels = ["(a) VORT_T","(b) OW_T","(c) LNAM"];

%Loop thru plots
for nTile = 1:numPanels
    ax = nexttile();

    %load the source of this panel
    thisSource = ["IO_" + PlotSources(nTile) + "_Comparisons.mat"];
    load(thisSource)

    %Grab maximum for line later
    [MaxVal,MaxInd] = max(mean(totalSims,1));

    %Plot individual lines
    plot(xToPlot,mean(totalNum,1),'--', 'linewidth', 6)
    hold on
    plot(xToPlot,mean(totalRatio,1),'--', 'linewidth', 6)
    plot(xToPlot,mean(totalSpatialC,1),'--', 'linewidth', 6)
    plot(xToPlot,mean(totalSpatialR,1), '--', 'linewidth', 6)
    plot(xToPlot,mean(totalSims,1)/100, 'linewidth', 9)
    plot([xToPlot(MaxInd) xToPlot(MaxInd)],[0 1], '--k', 'linewidth', 6)

    %labels and limits
    title(["Optimized at : " + num2str(xToPlot(MaxInd))])
    xlim([min(xToPlot) max(xToPlot)])
    ylim([0 1])
    ylabel("Scores")
    xlabel(xLabels(nTile))
    set(ax,'fontsize',44);
    grid on
    set(gca,'linew',5)
    if (nTile == numPanels)
        legend("Number Error","Ratio Error", "Spatial Correct Error", "Spatial Incorrect Error",...
            "Aggregate Similarity Score", "Maximum Score", 'Location', 'northeastoutside')
    end
    xTxT = min(xToPlot)+0.05*(max(xToPlot)-min(xToPlot));
    yTxT = 0.95;
    text(xTxT, yTxT, panelLabels(nTile), 'FontSize', 40)
    
end

%% Save and export
filenamestring = ([basepath + "/FIGURES/FINAL_RENDER/" + titlestring + ".tiff"]);
filename2save = char(filenamestring);
%we use export_fig because it automatically crops and is nice
export_fig(filename2save,'-m1.5','-a4','-opengl','-nocrop'); %saves to local directory as PNG
