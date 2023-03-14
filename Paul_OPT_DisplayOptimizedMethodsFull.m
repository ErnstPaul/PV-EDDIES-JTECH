%% Paul_OPT_DisplayOptimizedMethodsFull.m (version 2.2)
%Author: Paul Ernst
%Date Created: 6/1/2022
%Date of Last Update: 1/12/2023
%Update History:
%PE 1/12/2023 - Cleared in final revision
%PE 9/23/2022 - Cleaned up code for final release
%PE 6/13/2022 - Added additional panels for other metrics
%PE 6/1/2022 - Created
%--------------------------------------
%Purpose: Displays the optimized methods plot (Figure 5)
%Inputs: Plot sources from Comparison_Main
%Outputs: 1 panel for each method combination
%% Inputs
titlestring = "OptimizedMethod";
PlotSource = ["OPT_FINAL"];
xLabels = ["PV_I_S_O","OW_W_A","OW_T","VORT_W_A","VORT_T","LNAM"]; %list of methods to find center: SSH, PV_ISO, VORT, MV, OW, LNAM
yLabels = ["PV_I_S_O","OW_W_A","OW_T","VORT_W_A","VORT_T"]; %list of methods to find edge: SSH, PV_ISO, VORT, OW
beforestr = ["Err_N_u_m", "Err_R_a_t_i_o", "Err_P_o_s", "Err_N_e_g", "Aggregate Similarity Score"];
panelstr = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)"];

load('DefaultData_OPT.mat') %See Paul_OPT_CreateDefaultSettings.m
Paul_OPT_LoadDefaultSettings(); %See corresponding function
load("whiteorangecmap.mat")
%% Plots
figure('units', 'normalized', 'outerposition', [0 0.0652777777777778 0.64765625 0.917361111111111], 'Color', 'w'); %fullscreen, invisible figure

%Colorbar axes, your mileage might vary
cax1 = [0.26 0.9; 0.06 0.5; 0.7 1; 0.008 0.04; 50 75];

%load the source of this panel
thisSource = ["IO_" + PlotSource + "_Comparisons.mat"];
fields = {totalNum, totalRatio,  totalSpatialC, totalSpatialR, totalSims};
t = tiledlayout(3,2,'TileSpacing','tight');

for i = 1:5

    ax = nexttile();

    if i ~= 6

        field2use = fields{i};

        %Grab maximum for highlight later
        if i~=5
            [MaxVal,MaxInd] = min(mean(field2use,1));
            pcolToPlot = reshape(mean(field2use), [length(yLabels),length(xLabels)]); % x = smooth, y = ratio
            pcolToPlot(isnan(pcolToPlot)) = 100;
        elseif i == 5
            [MaxVal,MaxInd] = max(mean(field2use,1));
            pcolToPlot = reshape(mean(field2use), [length(yLabels),length(xLabels)]); % x = smooth, y = ratio
            pcolToPlot(isnan(pcolToPlot)) = 0;
        end
        [b, a] = ind2sub([length(yLabels),length(xLabels)],MaxInd);
        winnerY(i) = xLabels(a);
        winnerX(i) = yLabels(b);

        %Plot individual lines and the ticks of X's for maximums
        ss = pcolor(pcolToPlot);
        ss.FaceColor = 'texturemap';
        shading flat
        hold on
        colormap(CustomColormap)
        colorbar
        caxis(cax1(i,:))
        xTint = 5/6;
        yTint = 4/5;
        xticks(1+xTint/2:xTint:6)
        yticks(1+yTint/2:yTint:5)
        xlim([1 6])
        ylim([1 5])
        xticklabels(xLabels);
        yticklabels(yLabels);
        scatter(a*xTint+xTint*14/20,b*yTint+7/6*.5,2000,'kx','linewidth',3)

        %labels and limits
        title([panelstr(i) + " " + beforestr(i)])
        ylabel("Edge Method")
        xlabel("Center Method")
        set(gca,'fontsize',20);
        set(gca,'linew',8)
        messageToPrint(i,1) = [beforestr(i) + " Winner: " + winnerY(i) + " Centers; " + winnerX(i) + " Edges"];
    end
    if i == 5
        annotation(gcf,'textbox',[0.55003235620983 0.0837359098228664 0.405938693247344 0.221417069243157],...
            'String',[messageToPrint],...
            'FitBoxToText','off', 'FontSize', 28, 'LineWidth', 4);
    end
end

%% Save and export
filenamestring = ([basepath + "/FIGURES/FINAL_RENDER/" + titlestring + ".tiff"]);
filename2save = char(filenamestring);
%we use export_fig because it automatically crops and is nice
export_fig(filename2save,'-m1.5','-a4','-opengl','-nocrop'); %saves to local directory as PNG
